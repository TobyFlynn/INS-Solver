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

void xy_to_rs(const double x, const double y, double &r, double &s, const double *cellX, const double *cellY) {
  double l2 = (cellY[1] - cellY[2]) * (x - cellX[2]) + (cellX[2] - cellX[1]) * (y - cellY[2]);
  l2 = l2 / ((cellY[1] - cellY[2]) * (cellX[0] - cellX[2]) + (cellX[2] - cellX[1]) * (cellY[0] - cellY[2]));
  double l3 = (cellY[2] - cellY[0]) * (x - cellX[2]) + (cellX[0] - cellX[2]) * (y - cellY[2]);
  l3 = l3 / ((cellY[1] - cellY[2]) * (cellX[0] - cellX[2]) + (cellX[2] - cellX[1]) * (cellY[0] - cellY[2]));
  double l1 = 1.0 - l2 - l3;
  s = 2.0 * l1 - 1.0;
  r = 2.0 * l3 - 1.0;
}

void rs_to_xy(const double r, const double s, double &x, double &y, const double *cellX, const double *cellY) {
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

void eval_grad_at_pt(const double r, const double s, const double *modal,
                     double &dr, double &ds) {
  double a = -1.0;
  if(s != 1.0)
    a = 2.0 * (1.0 + r) / (1.0 - s) - 1.0;
  double b = s;
  arma::vec a_(1);
  arma::vec b_(1);
  a_[0] = a;
  b_[0] = b;

  dr = 0.0;
  ds = 0.0;
  int modal_ind = 0;
  arma::vec dr_v, ds_v;
  for(int x_ = 0; x_ <= 3; x_++) {
    for(int y_ = 0; y_ <= 3 - x_; y_++) {
      DGUtils::gradSimplex2DP(a_, b_, x_, y_, dr_v, ds_v);
      dr += modal[modal_ind] * dr_v[0];
      ds += modal[modal_ind++] * ds_v[0];
    }
  }
}

void eval_hessian_at_pt(const double r, const double s, const double *modal,
                        double &dr2, double &drs, double &ds2) {
  double a = -1.0;
  if(s != 1.0)
    a = 2.0 * (1.0 + r) / (1.0 - s) - 1.0;
  double b = s;
  arma::vec a_(1);
  arma::vec b_(1);
  a_[0] = a;
  b_[0] = b;

  dr2 = 0.0;
  drs = 0.0;
  ds2 = 0.0;
  int modal_ind = 0;
  arma::vec dr2_v, drs_v, ds2_v;
  for(int x_ = 0; x_ <= 3; x_++) {
    for(int y_ = 0; y_ <= 3 - x_; y_++) {
      DGUtils::hessianSimplex2DP(a_, b_, x_, y_, dr2_v, drs_v, ds2_v);
      dr2 += modal[modal_ind] * dr2_v[0];
      drs += modal[modal_ind] * drs_v[0];
      ds2 += modal[modal_ind++] * ds2_v[0];
    }
  }
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
                   const double *cellX, const double *cellY, const double rx, const double sx, const double ry, const double sy) {
  double lambda = 0.0;
  double node_r, node_s, pt_r, pt_s;
  xy_to_rs(node_x, node_y, node_r, node_s, cellX, cellY);
  xy_to_rs(closest_pt_x, closest_pt_y, pt_r, pt_s, cellX, cellY);
  double init_r = pt_r;
  double init_s = pt_s;
  double pt_x = closest_pt_x;
  double pt_y = closest_pt_y;
  double init_x = closest_pt_x;
  double init_y = closest_pt_y;
  const double h = 0.00288675;
  for(int step = 0; step < 10; step++) {
    double pt_r_old = pt_r;
    double pt_s_old = pt_s;
    double pt_x_old = pt_x;
    double pt_y_old = pt_y;
    // Evaluate surface and gradient at current guess
    double surface = eval_at_pt(pt_r, pt_s, s_modal);
    double surface_dr, surface_ds;
    eval_grad_at_pt(pt_r, pt_s, s_modal, surface_dr, surface_ds);
    const double surface_dx = rx * surface_dr + sx * surface_ds;
    const double surface_dy = ry * surface_dr + sy * surface_ds;
    // double surface_dr = eval_at_pt(pt_r, pt_s, dsdr_modal);
    // double surface_ds = eval_at_pt(pt_r, pt_s, dsds_modal);
    // Evaluate Hessian
    double hessian_tmp[3];
    eval_hessian_at_pt(pt_r, pt_s, s_modal, hessian_tmp[0], hessian_tmp[1], hessian_tmp[2]);
    double hessian[3];
    hessian[0] = rx * rx * hessian_tmp[0] + sx * sx * hessian_tmp[2];
    hessian[1] = rx * ry * hessian_tmp[1] + sx * sy * hessian_tmp[1];
    hessian[2] = ry * ry * hessian_tmp[0] + sy * sy * hessian_tmp[2];
    // hessian[0] = eval_at_pt(pt_r, pt_s, dsdr2_modal);
    // hessian[1] = eval_at_pt(pt_r, pt_s, dsdrs_modal);
    // hessian[2] = eval_at_pt(pt_r, pt_s, dsds2_modal);

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

    if(arma::cond(hessianf) > 1e3) {
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
    if((init_x - pt_x) * (init_x - pt_x) + (init_y - pt_y) * (init_y - pt_y) > h * h) {
      pt_x = pt_x_old;
      pt_y = pt_y_old;
      break;
    }

    // Converged, no more steps required
    if((pt_x_old - pt_x) * (pt_x_old - pt_x) + (pt_y_old - pt_y) * (pt_y_old - pt_y) < 1e-12)
      break;
  }

  // rs_to_xy(pt_r, pt_s, closest_pt_x, closest_pt_y, cellX, cellY);
  closest_pt_x = pt_x;
  closest_pt_y = pt_y;
}
