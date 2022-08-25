#include "ls.h"

#include "op_seq.h"

#include <limits>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <map>
#include <set>

#include "dg_constants.h"
#include "dg_utils.h"
#include "dg_blas_calls.h"
#include "dg_op2_blas.h"
#include "dg_operators.h"
#include "dg_compiler_defs.h"

#include "kd_tree_mpi_naive.h"
#include "timing.h"
#include "ls_reinit_poly.h"
#include "ls_reinit_poly_eval_cuda.h"
#include "utils.h"

#define THREADS_PER_BLOCK 256

extern Timing *timer;

__device__ void swap(double &a, double &b) {
  double c = a;
  a = b;
  b = c;
}

__device__ bool newtoncp_gepp(double A[][3], double *b) {
  for(int i = 0; i < 3; ++i) {
    int j = i;
    for(int k = i + 1; k < 3; ++k)
      if (fabs(A[k][i]) > fabs(A[j][i]))
        j = k;
    if (j != i) {
      for (int k = 0; k < 3; ++k)
        swap(A[i][k], A[j][k]);
      swap(b[i], b[j]);
    }

    if (fabs(A[i][i]) < 1.0e-8)
      return false;

    double fac = 1.0 / A[i][i];
    for (int j = i + 1; j < 3; ++j)
      A[j][i] *= fac;

    for (int j = i + 1; j < 3; ++j) {
      for (int k = i + 1; k < 3; ++k)
        A[j][k] -= A[j][i]*A[i][k];
      b[j] -= A[j][i]*b[i];
    }
  }

  for (int i = 3 - 1; i >= 0; --i) {
    double sum = 0.0;
    for (int j = i + 1; j < 3; ++j)
      sum += A[i][j]*b[j];
    b[i] = (b[i] - sum) / A[i][i];
  }

  return true;
}

__global__ void newton_kernel(const int numPts, const int *poly_ind, double *s,
                              const double *x, const double *y,
                              double *closest_x, double *closest_y,
                              PolyEval *pe, const double h) {
  // Assume only have threads in x dim
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if(tid >= numPts) return;

  const double x_l         = x[tid];
  const double y_l         = y[tid];
  const double closest_x_l = closest_x[tid];
  const double closest_y_l = closest_y[tid];

  double lambda = 0.0;
  double pt_x = closest_x_l;
  double pt_y = closest_y_l;
  double init_x = closest_x_l;
  double init_y = closest_y_l;

  for(int step = 0; step < 10; step++) {
    double pt_x_old = pt_x;
    double pt_y_old = pt_y;
    // Evaluate surface and gradient at current guess
    double surface = pe->val_at(poly_ind[tid], pt_x, pt_y);
    double surface_dx, surface_dy;
    pe->grad_at(poly_ind[tid], pt_x, pt_y, surface_dx, surface_dy);
    // Evaluate Hessian
    double hessian[3];
    pe->hessian_at(poly_ind[tid], pt_x, pt_y, hessian[0], hessian[1], hessian[2]);

    // Check if |nabla(surface)| = 0, if so then return
    double gradsqrnorm = surface_dx * surface_dx + surface_dy * surface_dy;
    if(gradsqrnorm < 1e-14)
      break;

    // Init lambda at first step
    if(step == 0)
      lambda = ((x_l - pt_x) * surface_dx + (y_l - pt_y) * surface_dy) / gradsqrnorm;

    // Gradient of functional
    double gradf[3];
    gradf[0] = pt_x - x_l + lambda * surface_dx;
    gradf[1] = pt_y - y_l + lambda * surface_dy;
    gradf[2] = surface;

    // Calculate Hessian of functional
    double hessianf[3][3];
    hessianf[0][0] = 1.0 + lambda * hessian[0];
    hessianf[0][1] = lambda * hessian[1]; hessianf[1][0] = hessianf[0][1];
    hessianf[0][2] = surface_dx; hessianf[2][0] = hessianf[0][2];

    hessianf[1][1] = 1.0 + lambda * hessian[2];
    hessianf[1][2] = surface_dy; hessianf[2][1] = hessianf[1][2];

    hessianf[2][2] = 0.0;

    if(!newtoncp_gepp(hessianf, gradf)) {
      double delta1_x = (surface / gradsqrnorm) * surface_dx;
      double delta1_y = (surface / gradsqrnorm) * surface_dy;
      lambda = ((x_l - pt_x) * surface_dx + (y_l - pt_y) * surface_dy) / gradsqrnorm;
      double delta2_x = pt_x - x_l + lambda * surface_dx;
      double delta2_y = pt_y - y_l + lambda * surface_dy;
      double msqr = delta2_x * delta2_x + delta2_y * delta2_y;
      if(msqr > 0.1 * h * 0.1 * h) {
        delta2_x *= 0.1 * h / sqrt(msqr);
        delta2_y *= 0.1 * h / sqrt(msqr);
      }
      pt_x -= delta1_x + delta2_x;
      pt_y -= delta1_y + delta2_y;
    } else {
      // Clamp update
      double msqr = gradf[0] * gradf[0] + gradf[1] * gradf[1];
      if(msqr > h * 0.5 * h * 0.5) {
        gradf[0] = gradf[0] * 0.5 * h / sqrt(msqr);
        gradf[1] = gradf[1] * 0.5 * h / sqrt(msqr);
        gradf[2] = gradf[2] * 0.5 * h / sqrt(msqr);
      }

      // Update guess
      pt_x -= gradf[0];
      pt_y -= gradf[1];
      lambda -= gradf[2];
    }

    // Gone outside the element, return
    // if((init_x - pt_x) * (init_x - pt_x) + (init_y - pt_y) * (init_y - pt_y) > h * h) {
    //   pt_x = pt_x_old;
    //   pt_y = pt_y_old;
    //   break;
    // }

    // Converged, no more steps required
    if((pt_x_old - pt_x) * (pt_x_old - pt_x) + (pt_y_old - pt_y) * (pt_y_old - pt_y) < 1e-18)
      break;
  }

  bool negative = s[tid] < 0.0;
  s[tid] = (pt_x - x_l) * (pt_x - x_l) + (pt_y - y_l) * (pt_y - y_l);
  s[tid] = sqrt(s[tid]);
  if(negative) s[tid] *= -1.0;
}

void newton_method(const int numPts, double *closest_x, double *closest_y, 
                   const double *x, const double *y, int *poly_ind, 
                   PolyEval *pe, double *s, const double h) {
  int threadsPerBlock, blocks;
  if(numPts < THREADS_PER_BLOCK) {
    threadsPerBlock = numPts;
    blocks = 1;
  } else {
    threadsPerBlock = THREADS_PER_BLOCK;
    blocks = numPts /  THREADS_PER_BLOCK;
    if(numPts % THREADS_PER_BLOCK != 0)
      blocks++;
  }

  newton_kernel<<<blocks, threadsPerBlock>>>(numPts, poly_ind, s, x, y,
                                             closest_x, closest_y, pe, h);
}

void LS::reinit_ls() {
  timer->startTimer("LS - Reinit");
  
  timer->startTimer("LS - Sample Interface");
  op_par_loop(sample_interface, "sample_interface", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(s,           -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(s_sample_x,  -1, OP_ID, LS_SAMPLE_NP, "double", OP_WRITE),
              op_arg_dat(s_sample_y,  -1, OP_ID, LS_SAMPLE_NP, "double", OP_WRITE));

  timer->endTimer("LS - Sample Interface");

  timer->startTimer("LS - Construct K-D Tree");
  const double *sample_pts_x = getOP2PtrHost(s_sample_x, OP_READ);
  const double *sample_pts_y = getOP2PtrHost(s_sample_y, OP_READ);

  KDTreeMPINaive kdtree(sample_pts_x, sample_pts_y, LS_SAMPLE_NP * mesh->cells->size, mesh, s);

  releaseOP2PtrHost(s_sample_x, OP_READ, sample_pts_x);
  releaseOP2PtrHost(s_sample_y, OP_READ, sample_pts_y);
  timer->endTimer("LS - Construct K-D Tree");

  const double *x_ptr = getOP2PtrHost(mesh->x, OP_READ);
  const double *y_ptr = getOP2PtrHost(mesh->y, OP_READ);

  double *closest_x, *closest_y;
  int *poly_ind;
  cudaMallocManaged(&closest_x, DG_NP * mesh->cells->size * sizeof(double));
  cudaMallocManaged(&closest_y, DG_NP * mesh->cells->size * sizeof(double));
  cudaMallocManaged(&poly_ind, DG_NP * mesh->cells->size * sizeof(int));

  timer->startTimer("LS - Query K-D Tree");
  #pragma omp parallel for
  for(int i = 0; i < DG_NP * mesh->cells->size; i++) {
    // Get closest sample point
    KDCoord tmp = kdtree.closest_point(x_ptr[i], y_ptr[i]);
    closest_x[i] = tmp.x;
    closest_y[i] = tmp.y;
    poly_ind[i]  = tmp.poly;
  }

  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);

  // Map of cell ind to polynomial approximations
  std::vector<PolyApprox> polys = kdtree.get_polys();
  PolyEval pe(polys);
  timer->endTimer("LS - Query K-D Tree");


  timer->startTimer("LS - Newton Method");
  const double *x_ptr_d = getOP2PtrDevice(mesh->x, OP_READ);
  const double *y_ptr_d = getOP2PtrDevice(mesh->y, OP_READ);
  double *s_ptr_d = getOP2PtrDevice(s, OP_RW);

  newton_method(DG_NP * mesh->cells->size, closest_x, closest_y, x_ptr_d,
                y_ptr_d, poly_ind, &pe, s_ptr_d, h);

  releaseOP2PtrDevice(s, OP_RW, s_ptr_d);
  releaseOP2PtrDevice(mesh->x, OP_READ, x_ptr_d);
  releaseOP2PtrDevice(mesh->y, OP_READ, y_ptr_d);

  cudaFree(closest_x);
  cudaFree(closest_y);
  cudaFree(poly_ind);
  timer->endTimer("LS - Newton Method");

  timer->endTimer("LS - Reinit");
}
