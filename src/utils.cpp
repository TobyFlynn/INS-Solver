#include "op_seq.h"

#include <memory>

#include "dg_utils.h"

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

void newton_kernel(const double node_x, const double node_y,
                   double &closest_pt_x, double &closest_pt_y,
                   const double *s_modal, const double *cellX,
                   const double *cellY, const double rx, const double sx,
                   const double ry, const double sy, const double h) {
  double lambda = 0.0;
  double node_r, node_s, pt_r, pt_s;
  DGUtils::global_xy_to_rs(node_x, node_y, node_r, node_s, cellX, cellY);
  DGUtils::global_xy_to_rs(closest_pt_x, closest_pt_y, pt_r, pt_s, cellX, cellY);
  double pt_x = closest_pt_x;
  double pt_y = closest_pt_y;
  double init_x = closest_pt_x;
  double init_y = closest_pt_y;
  for(int step = 0; step < 10; step++) {
    double pt_x_old = pt_x;
    double pt_y_old = pt_y;
    // Evaluate surface and gradient at current guess
    double surface = DGUtils::val_at_pt(pt_r, pt_s, s_modal);
    double surface_dr, surface_ds;
    DGUtils::grad_at_pt(pt_r, pt_s, s_modal, surface_dr, surface_ds);
    const double surface_dx = rx * surface_dr + sx * surface_ds;
    const double surface_dy = ry * surface_dr + sy * surface_ds;
    // Evaluate Hessian
    double hessian_tmp[3];
    DGUtils::hessian_at_pt(pt_r, pt_s, s_modal, hessian_tmp[0], hessian_tmp[1], hessian_tmp[2]);
    double hessian[3];
    hessian[0] = rx * rx * hessian_tmp[0] + sx * sx * hessian_tmp[2];
    hessian[1] = rx * ry * hessian_tmp[1] + sx * sy * hessian_tmp[1];
    hessian[2] = ry * ry * hessian_tmp[0] + sy * sy * hessian_tmp[2];

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

    if(true || arma::cond(hessianf) > 1e3) {
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
      arma::vec ans = arma::solve(hessianf, gradf);

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
    if((pt_x_old - pt_x) * (pt_x_old - pt_x) + (pt_y_old - pt_y) * (pt_y_old - pt_y) < 1e-12)
      break;

    DGUtils::global_xy_to_rs(pt_x, pt_y, pt_r, pt_s, cellX, cellY);
  }

  closest_pt_x = pt_x;
  closest_pt_y = pt_y;
}

void newton_method(const int numPts, double *s, double *closest_x,
                   double *closest_y, int *cell_ind, const double *x,
                   const double *y, const double *s_modal, const double *cell_x,
                   const double *cell_y, const double *rx, const double *sx,
                   const double *ry, const double *sy, const double h) {
  #pragma omp parallel for
  for(int i = 0; i < numPts; i++) {
    const int ind_np = cell_ind[i] * DG_NP;
    const int ind_3  = cell_ind[i] * 3;
    newton_kernel(x[i], y[i], closest_x[i], closest_y[i], &s_modal[ind_np],
                  &cell_x[ind_3], &cell_y[ind_3], rx[ind_np], sx[ind_np],
                  ry[ind_np], sy[ind_np], h);
    double s_ = s[i];
    bool negative = s[i] < 0.0;
    s[i] = (closest_x[i] - x[i]) * (closest_x[i] - x[i]) + (closest_y[i] - y[i]) * (closest_y[i] - y[i]);
    s[i] = sqrt(s[i]);
    if(negative) s[i] *= -1.0;
    if(fabs(s_ - s[i]) < 1e-5) s[i] = s_;
  }
}
