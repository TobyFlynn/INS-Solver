/*

#include "solvers/2d/ls_solver.h"

#include <vector>

#ifdef INS_MPI
#include "ls_utils/2d/kd_tree_mpi.h"
#else
#include "ls_utils/2d/kd_tree.h"
#endif
#include "op2_utils.h"
#include "timing.h"
#include "ls_utils/2d/ls_reinit_poly.h"
#include "ls_utils/2d/ls_reinit_poly_eval_device.h"

#define THREADS_PER_BLOCK 256

extern Timing *timer;

__device__ void swap(DG_FP &a, DG_FP &b) {
  DG_FP c = a;
  a = b;
  b = c;
}

__device__ bool newtoncp_gepp(DG_FP A[][3], DG_FP *b) {
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

    DG_FP fac = 1.0 / A[i][i];
    for (int j = i + 1; j < 3; ++j)
      A[j][i] *= fac;

    for (int j = i + 1; j < 3; ++j) {
      for (int k = i + 1; k < 3; ++k)
        A[j][k] -= A[j][i]*A[i][k];
      b[j] -= A[j][i]*b[i];
    }
  }

  for (int i = 3 - 1; i >= 0; --i) {
    DG_FP sum = 0.0;
    for (int j = i + 1; j < 3; ++j)
      sum += A[i][j]*b[j];
    b[i] = (b[i] - sum) / A[i][i];
  }

  return true;
}

__global__ void newton_kernel(const int numPts, const int *poly_ind, DG_FP *s,
                              const DG_FP *x, const DG_FP *y,
                              DG_FP *closest_x, DG_FP *closest_y,
                              PolyEval *pe, const DG_FP h, const int stride) {
  // Assume only have threads in x dim
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if(tid >= numPts) return;

  #ifdef DG_OP2_SOA
  const int dat_ind = tid;
  #else
  const int set_size = numPts / DG_NP;
  const int dat_ind = tid % set_size + stride * (tid / set_size);
  #endif

  const DG_FP x_l         = x[dat_ind];
  const DG_FP y_l         = y[dat_ind];
  const DG_FP closest_x_l = closest_x[tid];
  const DG_FP closest_y_l = closest_y[tid];

  DG_FP lambda = 0.0;
  DG_FP pt_x = closest_x_l;
  DG_FP pt_y = closest_y_l;
  DG_FP init_x = closest_x_l;
  DG_FP init_y = closest_y_l;

  for(int step = 0; step < 10; step++) {
    DG_FP pt_x_old = pt_x;
    DG_FP pt_y_old = pt_y;
    // Evaluate surface and gradient at current guess
    DG_FP surface = pe->val_at(poly_ind[tid], pt_x, pt_y);
    DG_FP surface_dx, surface_dy;
    pe->grad_at(poly_ind[tid], pt_x, pt_y, surface_dx, surface_dy);
    // Evaluate Hessian
    DG_FP hessian[3];
    pe->hessian_at(poly_ind[tid], pt_x, pt_y, hessian[0], hessian[1], hessian[2]);

    // Check if |nabla(surface)| = 0, if so then return
    DG_FP gradsqrnorm = surface_dx * surface_dx + surface_dy * surface_dy;
    if(gradsqrnorm < 1e-14)
      break;

    // Init lambda at first step
    if(step == 0)
      lambda = ((x_l - pt_x) * surface_dx + (y_l - pt_y) * surface_dy) / gradsqrnorm;

    // Gradient of functional
    DG_FP gradf[3];
    gradf[0] = pt_x - x_l + lambda * surface_dx;
    gradf[1] = pt_y - y_l + lambda * surface_dy;
    gradf[2] = surface;

    // Calculate Hessian of functional
    DG_FP hessianf[3][3];
    hessianf[0][0] = 1.0 + lambda * hessian[0];
    hessianf[0][1] = lambda * hessian[1]; hessianf[1][0] = hessianf[0][1];
    hessianf[0][2] = surface_dx; hessianf[2][0] = hessianf[0][2];

    hessianf[1][1] = 1.0 + lambda * hessian[2];
    hessianf[1][2] = surface_dy; hessianf[2][1] = hessianf[1][2];

    hessianf[2][2] = 0.0;

    if(!newtoncp_gepp(hessianf, gradf)) {
      DG_FP delta1_x = (surface / gradsqrnorm) * surface_dx;
      DG_FP delta1_y = (surface / gradsqrnorm) * surface_dy;
      lambda = ((x_l - pt_x) * surface_dx + (y_l - pt_y) * surface_dy) / gradsqrnorm;
      DG_FP delta2_x = pt_x - x_l + lambda * surface_dx;
      DG_FP delta2_y = pt_y - y_l + lambda * surface_dy;
      DG_FP msqr = delta2_x * delta2_x + delta2_y * delta2_y;
      if(msqr > 0.1 * h * 0.1 * h) {
        delta2_x *= 0.1 * h / sqrt(msqr);
        delta2_y *= 0.1 * h / sqrt(msqr);
      }
      pt_x -= delta1_x + delta2_x;
      pt_y -= delta1_y + delta2_y;
    } else {
      // Clamp update
      DG_FP msqr = gradf[0] * gradf[0] + gradf[1] * gradf[1];
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

  bool negative = s[dat_ind] < 0.0;
  s[dat_ind] = (pt_x - x_l) * (pt_x - x_l) + (pt_y - y_l) * (pt_y - y_l);
  s[dat_ind] = sqrt(s[dat_ind]);
  if(negative) s[dat_ind] *= -1.0;
}

void newton_method(const int numPts, DG_FP *closest_x, DG_FP *closest_y,
                   const DG_FP *x, const DG_FP *y, int *poly_ind,
                   PolyEval *pe, DG_FP *s, const DG_FP h) {
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

  PolyEval *pe_d;
  cudaMallocManaged(&pe_d, sizeof(PolyEval));
  memcpy(pe_d, pe, sizeof(PolyEval));

  newton_kernel<<<blocks, threadsPerBlock>>>(numPts, poly_ind, s, x, y,
                                             closest_x, closest_y, pe_d, h);
  cudaFree(pe_d);
}

void LevelSetSolver2D::reinitLS() {
  timer->startTimer("LS - Reinit");

  sampleInterface();

  timer->startTimer("LS - Construct K-D Tree");
  const DG_FP *sample_pts_x = getOP2PtrHost(s_sample_x, OP_READ);
  const DG_FP *sample_pts_y = getOP2PtrHost(s_sample_y, OP_READ);

  #ifdef INS_MPI
  KDTreeMPI kdtree(sample_pts_x, sample_pts_y, LS_SAMPLE_NP * mesh->cells->size, mesh, s);
  #else
  KDTree kdtree(sample_pts_x, sample_pts_y, LS_SAMPLE_NP * mesh->cells->size, mesh, s);
  #endif

  releaseOP2PtrHost(s_sample_x, OP_READ, sample_pts_x);
  releaseOP2PtrHost(s_sample_y, OP_READ, sample_pts_y);
  timer->endTimer("LS - Construct K-D Tree");

  const DG_FP *x_ptr = getOP2PtrHost(mesh->x, OP_READ);
  const DG_FP *y_ptr = getOP2PtrHost(mesh->y, OP_READ);

  DG_FP *closest_x, *closest_y;
  int *poly_ind;
  cudaMallocManaged(&closest_x, DG_NP * mesh->cells->size * sizeof(DG_FP));
  cudaMallocManaged(&closest_y, DG_NP * mesh->cells->size * sizeof(DG_FP));
  cudaMallocManaged(&poly_ind, DG_NP * mesh->cells->size * sizeof(int));

  timer->startTimer("LS - Query K-D Tree");

  #ifdef INS_MPI

  const DG_FP *s_ptr = getOP2PtrHost(s, OP_READ);
  int num_pts_to_reinit = 0;
  std::vector<DG_FP> x_vec, y_vec;
  for(int i = 0; i < DG_NP * mesh->cells->size; i++) {
    int start_ind = (i / DG_NP) * DG_NP;
    bool reinit = false;
    for(int j = 0; j < DG_NP; j++) {
      if(fabs(s_ptr[start_ind + j]) < 0.05) {
        reinit = true;
      }
    }
    if(reinit) {
      num_pts_to_reinit++;
      x_vec.push_back(x_ptr[i]);
      y_vec.push_back(y_ptr[i]);
    }
  }

  std::vector<DG_FP> cx_vec(num_pts_to_reinit), cy_vec(num_pts_to_reinit);
  std::vector<int> p_vec(num_pts_to_reinit);

  kdtree.closest_point(num_pts_to_reinit, x_vec.data(), y_vec.data(), cx_vec.data(), cy_vec.data(), p_vec.data());

  int count = 0;
  for(int i = 0; i < DG_NP * mesh->cells->size; i++) {
    int start_ind = (i / DG_NP) * DG_NP;
    bool reinit = false;
    for(int j = 0; j < DG_NP; j++) {
      if(fabs(s_ptr[start_ind + j]) < 0.05) {
        reinit = true;
      }
    }
    if(reinit) {
      closest_x[i] = cx_vec[count];
      closest_y[i] = cy_vec[count];
      poly_ind[i] = p_vec[count];
      count++;
    }
  }

  releaseOP2PtrHost(s, OP_READ, s_ptr);

  #else

  #pragma omp parallel for
  for(int i = 0; i < DG_NP * mesh->cells->size; i++) {
    // Get closest sample point
    KDCoord tmp = kdtree.closest_point(x_ptr[i], y_ptr[i]);
    closest_x[i] = tmp.x;
    closest_y[i] = tmp.y;
    // Convert OP2 cell ind to PolyEval ind (minimise indirection for CUDA)
    poly_ind[i]  = tmp.poly;
  }

  #endif

  timer->endTimer("LS - Query K-D Tree");

  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);

  // Map of cell ind to polynomial approximations
  timer->startTimer("LS - Construct Poly Eval");
  std::vector<PolyApprox> polys = kdtree.get_polys();
  PolyEval pe(polys);
  timer->endTimer("LS - Construct Poly Eval");

  #ifdef DG_OP2_SOA
  throw std::runtime_error("2D LS Reinit not implemented for SoA");
  #endif

  timer->startTimer("LS - Newton Method");
  const DG_FP *x_ptr_d = getOP2PtrDevice(mesh->x, OP_READ);
  const DG_FP *y_ptr_d = getOP2PtrDevice(mesh->y, OP_READ);
  DG_FP *s_ptr_d = getOP2PtrDevice(s, OP_RW);

  op_arg args[] = {
    op_arg_dat(mesh->x, -1, OP_ID, dat->dim, DG_FP_STR, acc)
  };
  const int stride_size = getSetSizeFromOpArg(&args[0]);

  // newton_method(DG_NP * mesh->cells->size, closest_x, closest_y, x_ptr_d,
  //               y_ptr_d, poly_ind, &pe, s_ptr_d, h, stride_size);

  releaseOP2PtrDevice(s, OP_RW, s_ptr_d);
  releaseOP2PtrDevice(mesh->x, OP_READ, x_ptr_d);
  releaseOP2PtrDevice(mesh->y, OP_READ, y_ptr_d);

  cudaFree(closest_x);
  cudaFree(closest_y);
  cudaFree(poly_ind);
  timer->endTimer("LS - Newton Method");
  timer->endTimer("LS - Reinit");
}
*/

#include "ls_reinit.cpp"
