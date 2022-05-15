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

double xy_to_rs(const double x, const double y, double &r, double &s, const double *cellX, const double *cellY) {
  s = (cellY[1] - cellY[0]) * (2.0 * x - cellX[1] - cellX[2]) - (cellX[1] - cellX[0]) * (2.0 * y - cellY[1] - cellY[2]);
  s = s / ((cellX[0] - cellX[1]) * (cellY[2] - cellY[0]) + (cellY[1] - cellY[0]) * (cellX[2] - cellX[0]));
  r = 2.0 * x - s * (cellX[2] - cellX[0]) - cellX[1] - cellX[2];
  r = r / (cellX[1] - cellX[0]);
}

double rs_to_xy(const double r, const double s, double &x, double &y, const double *cellX, const double *cellY) {
  x = 0.5 * (-(r + s) * cellX[0] + (1.0 + r) * cellX[1] + (1.0 + s) * cellX[2]);
  y = 0.5 * (-(r + s) * cellY[0] + (1.0 + r) * cellY[1] + (1.0 + s) * cellY[2]);
}

double eval_at_pt(const double r, const double s, const double *modal) {
  double a = -1.0;
  if(s != 1.0)
    a = 2.0 * (1.0 + r) / (1.0 - s) - 1.0;
  double b = s;
  arma::vec a_(1);
  arma::vec b_(1);
  a_[0] = a;
  b_[0] = b;

  double new_val = 0.0;
  int modal_ind = 0;
  for(int x_ = 0; x_ <= 3; x_++) {
    for(int y_ = 0; y_ <= 3 - x_; y_++) {
      arma::vec res = DGUtils::simplex2DP(a_, b_, x_, y_);
      new_val += modal[modal_ind++] * res[0];
    }
  }

  return new_val;
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

void newton_method(const double node_x, const double node_y,
                   double &closest_pt_x, double &closest_pt_y,
                   const double *s_modal, const double *dsdr_modal,
                   const double *dsds_modal, const double *dsdr2_modal,
                   const double *dsdrs_modal, const double *dsds2_modal,
                   const double *cellX, const double *cellY) {
  double lambda = 0.0;
  double node_r, node_s, pt_r, pt_s;
  xy_to_rs(node_x, node_y, node_r, node_s, cellX, cellY);
  xy_to_rs(closest_pt_x, closest_pt_y, pt_r, pt_s, cellX, cellY);
  double init_r = pt_r;
  double init_s = pt_s;
  for(int step = 0; step < 10; step++) {
    double pt_r_old = pt_r;
    double pt_s_old = pt_s;
    // Evaluate surface and gradient at current guess
    double surface    = eval_at_pt(pt_r, pt_s, s_modal);
    double surface_dr = eval_at_pt(pt_r, pt_s, dsdr_modal);
    double surface_ds = eval_at_pt(pt_r, pt_s, dsds_modal);
    // Evaluate Hessian
    double hessian[3];
    hessian[0] = eval_at_pt(pt_r, pt_s, dsdr2_modal);
    hessian[1] = eval_at_pt(pt_r, pt_s, dsdrs_modal);
    hessian[2] = eval_at_pt(pt_r, pt_s, dsds2_modal);

    // Check if |nabla(surface)| = 0, if so then return
    double gradsqrnorm = surface_dr * surface_dr + surface_ds * surface_ds;
    if(gradsqrnorm < 1e-10)
      return;

    // Init lambda at first step
    if(step == 0)
      lambda = ((node_r - pt_r) * surface_dr + (node_s - pt_s) * surface_ds) / gradsqrnorm;

    // Gradient of functional
    arma::vec gradf(3);
    gradf(0) = pt_r - node_r + lambda * surface_dr;
    gradf(1) = pt_s - node_s + lambda * surface_ds;
    gradf(2) = surface;

    // Calculate Hessian of functional
    arma::mat hessianf(3, 3);
    hessianf(0, 0) = 1.0 + lambda * hessian[0];
    hessianf(0, 1) = lambda * hessian[1]; hessianf(1, 0) = hessianf(0, 1);
    hessianf(0, 2) = surface_dr; hessianf(2, 0) = hessianf(0, 2);

    hessianf(1, 1) = 1.0 + lambda * hessian[2];
    hessianf(1, 2) = surface_ds; hessianf(2, 1) = hessianf(1, 2);

    hessianf(2, 2) = 0.0;

    if(arma::cond(hessianf) > 1e5)
      break;

    arma::vec ans = arma::solve(hessianf, gradf);

    // Clamp update
    double msqr = ans(0) * ans(0) + ans(1) * ans(1);
    if(msqr > 1.5 * 0.5 * 1.5 * 0.5)
      ans = ans * 0.5 * 1.5 / (msqr * msqr);

    // Update guess
    pt_r -= ans(0);
    pt_s -= ans(1);
    lambda -= ans(2);

    // Gone outside the element, return
    if((init_r - pt_r) * (init_r - pt_r) + (init_s - pt_s) * (init_s - pt_s) > 1.5 * 1.5)
      return;

    // Converged, no more steps required
    if((pt_r_old - pt_r) * (pt_r_old - pt_r) + (pt_s_old - pt_s) * (pt_s_old - pt_s) < 1e-8)
      break;
  }

  rs_to_xy(pt_r, pt_s, closest_pt_x, closest_pt_y, cellX, cellY);
}
