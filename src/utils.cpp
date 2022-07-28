#include "op_seq.h"

#include <memory>
#include <iostream>

#include "dg_utils.h"

#include "ls_reinit_stencil.h"

double *getOP2Array(op_dat dat) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_SUB_CELLS, "double", OP_READ)
  };
  op_mpi_halo_exchanges(dat->set, 1, args);
  double *res = (double *)malloc(dat->set->size * DG_SUB_CELLS * sizeof(double));
  memcpy(res, dat->data, dat->set->size * DG_SUB_CELLS * sizeof(double));
  op_mpi_set_dirtybit(1, args);
  return res;
}

bool is_point_in_cell(const double x, const double y, const double *cellX, const double *cellY) {
  double ABx = cellX[1] - cellX[0];
  double ABy = cellY[1] - cellY[0];

  double APx = x - cellX[0];
  double APy = y - cellY[0];

  double BCx = cellX[2] - cellX[1];
  double BCy = cellY[2] - cellY[1];

  double BPx = x - cellX[1];
  double BPy = y - cellY[1];

  double CAx = cellX[0] - cellX[2];
  double CAy = cellY[0] - cellY[2];

  double CPx = x - cellX[2];
  double CPy = y - cellY[2];

  double AB_AP = ABx * APy - ABy * APx;
  double BC_BP = BCx * BPy - BCy * BPx;
  double CA_CP = CAx * CPy - CAy * CPx;

  bool zero0 = AB_AP == 0.0 || AB_AP == -0.0;
  bool zero1 = BC_BP == 0.0 || BC_BP == -0.0;
  bool zero2 = CA_CP == 0.0 || CA_CP == -0.0;

  if(zero0 || zero1 || zero2) {
    // One zero means on an edge of the triangle
    // Two zeros means on a vertex of the triangle
    if(zero0 && !zero1 && !zero2) {
      return (BC_BP > 0.0 && CA_CP > 0.0) || (BC_BP < 0.0 && CA_CP < 0.0);
    } else if(!zero0 && zero1 && !zero2) {
      return (AB_AP > 0.0 && CA_CP > 0.0) || (AB_AP < 0.0 && CA_CP < 0.0);
    } else if(!zero0 && !zero1 && zero2) {
      return (AB_AP > 0.0 && BC_BP > 0.0) || (AB_AP < 0.0 && BC_BP < 0.0);
    } else if(zero0 && zero1 && !zero2) {
      return true;
    } else if(zero0 && !zero1 && zero2) {
      return true;
    } else if(!zero0 && zero1 && zero2) {
      return true;
    } else {
      return false;
    }
  }

  if(AB_AP > 0.0 && BC_BP > 0.0 && CA_CP > 0.0) {
    return true;
  } else if(AB_AP < 0.0 && BC_BP < 0.0 && CA_CP < 0.0) {
    return true;
  } else {
    return false;
  }
}

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
                   const int cell_ind, op_map edge_map, op_dat x_dat, op_dat y_dat,
                   op_dat s_dat, const double h) {
  double lambda = 0.0;
  bool converged = false;
  double pt_x = closest_pt_x;
  double pt_y = closest_pt_y;
  double init_x = closest_pt_x;
  double init_y = closest_pt_y;

  PolyApprox p(2, cell_ind, edge_map, x_dat, y_dat, s_dat);

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
                   int *cell_ind, op_map edge_map, op_dat x_dat, op_dat y_dat, op_dat s_dat, const double h) {
  int numNonConv = 0;
  int numReinit = 0;
  double *s_ptr = (double *)s_dat->data;
  double *s_new = (double *)calloc(numPts, sizeof(double));
  #pragma omp parallel for
  for(int i = 0; i < numPts; i++) {
    int start_ind = (i / DG_NP) * DG_NP;
    bool reinit = false;
    for(int j = 0; j < DG_NP; j++) {
      if(fabs(s_ptr[start_ind + j]) < 0.05) {
        reinit = true;
      }
    }

    s_new[i] = s_ptr[i];

    if(reinit) {
      bool tmp = newton_kernel(closest_x[i], closest_y[i], x[i], y[i], cell_ind[i], 
                               edge_map, x_dat, y_dat, s_dat, h);
      if(tmp) {
        bool negative = s_ptr[i] < 0.0;
        s_new[i] = (closest_x[i] - x[i]) * (closest_x[i] - x[i]) + (closest_y[i] - y[i]) * (closest_y[i] - y[i]);
        s_new[i] = sqrt(s_new[i]);
        if(negative) s_new[i] *= -1.0;
      }
      if(!tmp) {
        #pragma omp atomic
        numNonConv++;
      }
      #pragma omp atomic
      numReinit++;
    }
  }

  memcpy(s_ptr, s_new, numPts * sizeof(double));

  free(s_new);
  
  if(numNonConv != 0 || numReinit == 0)
    std::cout << numNonConv << " non-converged points out of " << numReinit << " points reinitialised" << std::endl;
}
