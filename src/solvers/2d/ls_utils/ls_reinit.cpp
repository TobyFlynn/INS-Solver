#include "solvers/2d/ls_solver.h"

#include <vector>
#include <memory>

#ifdef INS_MPI
#include "ls_utils/2d/kd_tree_mpi.h"
#else
#include "ls_utils/2d/kd_tree.h"
#endif
#include "timing.h"
#include "ls_utils/2d/ls_reinit_poly.h"
#include "op2_utils.h"

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

    if (std::abs(A(i,i)) < 1.0e4*std::numeric_limits<DG_FP>::epsilon())
      return false;

    DG_FP fac = 1.0 / A(i,i);
    for (int j = i + 1; j < 3; ++j)
      A(j,i) *= fac;

    for (int j = i + 1; j < 3; ++j) {
      for (int k = i + 1; k < 3; ++k)
        A(j,k) -= A(j,i)*A(i,k);
      b(j) -= A(j,i)*b(i);
    }
  }

  for (int i = 3 - 1; i >= 0; --i) {
    DG_FP sum = 0.0;
    for (int j = i + 1; j < 3; ++j)
      sum += A(i,j)*b(j);
    b(i) = (b(i) - sum) / A(i,i);
  }

  return true;
}

bool newton_kernel(DG_FP &closest_pt_x, DG_FP &closest_pt_y,
                   const DG_FP node_x, const DG_FP node_y,
                   PolyApprox &p, const DG_FP h) {
  DG_FP lambda = 0.0;
  bool converged = false;
  DG_FP pt_x = closest_pt_x;
  DG_FP pt_y = closest_pt_y;
  DG_FP init_x = closest_pt_x;
  DG_FP init_y = closest_pt_y;

  for(int step = 0; step < 10; step++) {
    DG_FP pt_x_old = pt_x;
    DG_FP pt_y_old = pt_y;
    // Evaluate surface and gradient at current guess
    DG_FP surface = p.val_at(pt_x, pt_y);
    DG_FP surface_dx, surface_dy;
    p.grad_at(pt_x, pt_y, surface_dx, surface_dy);
    // Evaluate Hessian
    DG_FP hessian[3];
    p.hessian_at(pt_x, pt_y, hessian[0], hessian[1], hessian[2]);

    // Check if |nabla(surface)| = 0, if so then return
    DG_FP gradsqrnorm = surface_dx * surface_dx + surface_dy * surface_dy;
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
      DG_FP delta1_x = (surface / gradsqrnorm) * surface_dx;
      DG_FP delta1_y = (surface / gradsqrnorm) * surface_dy;
      lambda = ((node_x - pt_x) * surface_dx + (node_y - pt_y) * surface_dy) / gradsqrnorm;
      DG_FP delta2_x = pt_x - node_x + lambda * surface_dx;
      DG_FP delta2_y = pt_y - node_y + lambda * surface_dy;
      DG_FP msqr = delta2_x * delta2_x + delta2_y * delta2_y;
      if(msqr > 0.1 * h * 0.1 * h) {
        delta2_x *= 0.1 * h / sqrt(msqr);
        delta2_y *= 0.1 * h / sqrt(msqr);
      }
      pt_x -= delta1_x + delta2_x;
      pt_y -= delta1_y + delta2_y;
    } else {
      arma::vec ans = gradf;

      // Clamp update
      DG_FP msqr = ans(0) * ans(0) + ans(1) * ans(1);
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

void newton_method(const int numPts, DG_FP *closest_x, DG_FP *closest_y,
                   const DG_FP *x, const DG_FP *y, int *poly_ind,
                   std::vector<PolyApprox> &polys, DG_FP *s, const DG_FP h,
                   const DG_FP reinit_width) {
  int numNonConv = 0;
  int numReinit = 0;

  #pragma omp parallel for
  for(int i = 0; i < numPts; i++) {
    int start_ind = (i / DG_NP) * DG_NP;
    bool reinit = false;
    for(int j = 0; j < DG_NP; j++) {
      if(fabs(s[start_ind + j]) < reinit_width) {
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

  double percent_non_converge = numReinit == 0 ? 0.0 : (double)numNonConv / (double)numReinit;
  if(percent_non_converge > 0.1)
    std::cout << percent_non_converge * 100.0 << "\% reinitialisation points did not converge" << std::endl;
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

  DG_FP *closest_x = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  DG_FP *closest_y = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  int *poly_ind     = (int *)calloc(DG_NP * mesh->cells->size, sizeof(int));

  timer->startTimer("LS - Query K-D Tree");
  #ifdef INS_MPI

  const DG_FP *s_ptr = getOP2PtrHost(s, OP_READ);
  int num_pts_to_reinit = 0;
  std::vector<DG_FP> x_vec, y_vec;
  for(int i = 0; i < DG_NP * mesh->cells->size; i++) {
    int start_ind = (i / DG_NP) * DG_NP;
    bool reinit = false;
    for(int j = 0; j < DG_NP; j++) {
      if(fabs(s_ptr[start_ind + j]) < reinit_width) {
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
      if(fabs(s_ptr[start_ind + j]) < reinit_width) {
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
    poly_ind[i]  = tmp.poly;
  }

  #endif
  timer->endTimer("LS - Query K-D Tree");

  // Map of cell ind to polynomial approximations
  std::vector<PolyApprox> polys = kdtree.get_polys();

  timer->startTimer("LS - Newton Method");
  DG_FP *surface_ptr = getOP2PtrHost(s, OP_RW);

  newton_method(DG_NP * mesh->cells->size, closest_x, closest_y, x_ptr, y_ptr,
                poly_ind, polys, surface_ptr, h, reinit_width);

  releaseOP2PtrHost(s, OP_RW, surface_ptr);

  free(closest_x);
  free(closest_y);
  free(poly_ind);

  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);
  timer->endTimer("LS - Newton Method");
  timer->endTimer("LS - Reinit");
}
