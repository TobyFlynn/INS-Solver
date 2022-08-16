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
#include "utils.h"

extern Timing *timer;

bool newtoncp_gepp(arma::mat &A, arma::vec &b) {
  for (int i = 0; i < 3; ++i) {
    int j = i;
    for (int k = i + 1; k < 3; ++k)
      if (std::abs(A(k,i)) > std::abs(A(j,i)))
        j = k;
    if (j != i) {
      for (int k = 0; k < 3; ++k)
        std::swap(A(i,k), A(j,k));
      std::swap(b(i), b(j));
    }

    if (std::abs(A(i,i)) < 1.0e4*std::numeric_limits<double>::epsilon())
      return false;

    double fac = 1.0 / A(i,i);
    for (int j = i + 1; j < 3; ++j)
      A(j,i) *= fac;

    for (int j = i + 1; j < 3; ++j) {
      for (int k = i + 1; k < 3; ++k)
        A(j,k) -= A(j,i)*A(i,k);
      b(j) -= A(j,i)*b(i);
    }
  }

  for (int i = 3 - 1; i >= 0; --i) {
    double sum = 0.0;
    for (int j = i + 1; j < 3; ++j)
      sum += A(i,j)*b(j);
    b(i) = (b(i) - sum) / A(i,i);
  }

  return true;
}

bool newton_kernel(double &closest_pt_x, double &closest_pt_y,
                   const double node_x, const double node_y,
                   PolyApprox &p, const double h) {
  double lambda = 0.0;
  bool converged = false;
  double pt_x = closest_pt_x;
  double pt_y = closest_pt_y;
  double init_x = closest_pt_x;
  double init_y = closest_pt_y;

  for(int step = 0; step < 10; step++) {
    double pt_x_old = pt_x;
    double pt_y_old = pt_y;
    // Evaluate surface and gradient at current guess
    double surface = p.val_at(pt_x, pt_y);
    double surface_dx, surface_dy;
    p.grad_at(pt_x, pt_y, surface_dx, surface_dy);
    // Evaluate Hessian
    double hessian[3];
    p.hessian_at(pt_x, pt_y, hessian[0], hessian[1], hessian[2]);
    
    // Check if |nabla(surface)| = 0, if so then return
    double gradsqrnorm = surface_dx * surface_dx + surface_dy * surface_dy;
    if(gradsqrnorm < 1e-14)
      break;

    // Init lambda at first step
    if(step == 0)
      lambda = ((node_x - pt_x) * surface_dx + (node_y - pt_y) * surface_dy) / gradsqrnorm;

    // Gradient of functional
    arma::vec gradf(3);
    gradf(0) = pt_x - node_x + lambda * surface_dx;
    gradf(1) = pt_y - node_y + lambda * surface_dy;
    gradf(2) = surface;

    // Calculate Hessian of functional
    arma::mat hessianf(3, 3);
    hessianf(0, 0) = 1.0 + lambda * hessian[0];
    hessianf(0, 1) = lambda * hessian[1]; hessianf(1, 0) = hessianf(0, 1);
    hessianf(0, 2) = surface_dx; hessianf(2, 0) = hessianf(0, 2);

    hessianf(1, 1) = 1.0 + lambda * hessian[2];
    hessianf(1, 2) = surface_dy; hessianf(2, 1) = hessianf(1, 2);

    hessianf(2, 2) = 0.0;

    if(!newtoncp_gepp(hessianf, gradf)) {
      double delta1_x = (surface / gradsqrnorm) * surface_dx;
      double delta1_y = (surface / gradsqrnorm) * surface_dy;
      lambda = ((node_x - pt_x) * surface_dx + (node_y - pt_y) * surface_dy) / gradsqrnorm;
      double delta2_x = pt_x - node_x + lambda * surface_dx;
      double delta2_y = pt_y - node_y + lambda * surface_dy;
      double msqr = delta2_x * delta2_x + delta2_y * delta2_y;
      if(msqr > 0.1 * h * 0.1 * h) {
        delta2_x *= 0.1 * h / sqrt(msqr);
        delta2_y *= 0.1 * h / sqrt(msqr);
      }
      pt_x -= delta1_x + delta2_x;
      pt_y -= delta1_y + delta2_y;
    } else {
      arma::vec ans = gradf;

      // Clamp update
      double msqr = ans(0) * ans(0) + ans(1) * ans(1);
      if(msqr > h * 0.5 * h * 0.5)
        ans = ans * 0.5 * h / sqrt(msqr);

      // Update guess
      pt_x -= ans(0);
      pt_y -= ans(1);
      lambda -= ans(2);
    }

    // Gone outside the element, return
    // if((init_x - pt_x) * (init_x - pt_x) + (init_y - pt_y) * (init_y - pt_y) > h * h) {
    //   pt_x = pt_x_old;
    //   pt_y = pt_y_old;
    //   break;
    // }

    // Converged, no more steps required
    if((pt_x_old - pt_x) * (pt_x_old - pt_x) + (pt_y_old - pt_y) * (pt_y_old - pt_y) < 1e-18) {
      converged = true;
      break;
    }
  }

  closest_pt_x = pt_x;
  closest_pt_y = pt_y;

  return converged;
}

void newton_method(const int numPts, double *closest_x, double *closest_y, const double *x, const double *y,
                   int *poly_ind, std::vector<PolyApprox> &polys, double *s, const double h) {
  int numNonConv = 0;
  int numReinit = 0;

  #pragma omp parallel for
  for(int i = 0; i < numPts; i++) {
    int start_ind = (i / DG_NP) * DG_NP;
    bool reinit = false;
    for(int j = 0; j < DG_NP; j++) {
      if(fabs(s[start_ind + j]) < 0.05) {
        reinit = true;
      }
    }

    if(reinit) {
      bool tmp = newton_kernel(closest_x[i], closest_y[i], x[i], y[i], polys[poly_ind[i]], h);
      if(tmp) {
        bool negative = s[i] < 0.0;
        s[i] = (closest_x[i] - x[i]) * (closest_x[i] - x[i]) + (closest_y[i] - y[i]) * (closest_y[i] - y[i]);
        s[i] = sqrt(s[i]);
        if(negative) s[i] *= -1.0;
      }
      if(!tmp) {
        #pragma omp atomic
        numNonConv++;
      }
      #pragma omp atomic
      numReinit++;
    }
  }
  
  if(numNonConv != 0)
    std::cout << numNonConv << " non-converged points out of " << numReinit << " points reinitialised" << std::endl;
}

void LS::reinit_ls() {
  timer->startTimer("LS - Reinit");
  timer->startTimer("LS - Sample Interface");
  op_par_loop(sample_interface, "sample_interface", mesh->cells,
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

  double *closest_x = (double *)calloc(DG_NP * mesh->cells->size, sizeof(double));
  double *closest_y = (double *)calloc(DG_NP * mesh->cells->size, sizeof(double));
  int *poly_ind     = (int *)calloc(DG_NP * mesh->cells->size, sizeof(int));

  timer->startTimer("LS - Query K-D Tree");
  #pragma omp parallel for
  for(int i = 0; i < DG_NP * mesh->cells->size; i++) {
    // Get closest sample point
    KDCoord tmp = kdtree.closest_point(x_ptr[i], y_ptr[i]);
    closest_x[i] = tmp.x;
    closest_y[i] = tmp.y;
    poly_ind[i]  = tmp.poly;
  }

  // Map of cell ind to polynomial approximations
  std::vector<PolyApprox> polys = kdtree.get_polys();
  timer->endTimer("LS - Query K-D Tree");


  timer->startTimer("LS - Newton Method");
  double *surface_ptr = getOP2PtrHost(s, OP_RW);

  newton_method(DG_NP * mesh->cells->size, closest_x, closest_y, x_ptr, y_ptr,
                poly_ind, polys, surface_ptr, h);

  releaseOP2PtrHost(s, OP_RW, surface_ptr);

  free(closest_x);
  free(closest_y);
  free(poly_ind);

  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);
  timer->endTimer("LS - Newton Method");
  timer->endTimer("LS - Reinit");
}
