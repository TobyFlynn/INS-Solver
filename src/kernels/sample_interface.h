#include "dg_utils.h"

inline double eval_at_pt(const double r, const double s, const double *modal) {
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

inline bool is_point_in_cell(const double x, const double y, const double *cellX, const double *cellY) {
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

inline void sample_interface(const double *h, const double *x, const double *y,
                             const double *surface, const double *s_modal,
                             const double *dsdx_modal, const double *dsdy_modal,
                             double *sample_x, double *sample_y) {
  // Check that the cell contains the interface
  bool positive0 = surface[0] > 0.0;
  bool interface = false;
  for(int i = 1; i < DG_NP; i++) {
    if(positive0 != surface[i] > 0.0)
      interface = true;
  }
  if(!interface)
    return;
  // Initial positions of sample points (in r-s coords)
  sample_x[0] = -0.55158350755530561;
  sample_y[0] = -0.55158350755530561;
  sample_x[1] = 0.10316701511061113;
  sample_y[1] = -0.55158350755530561;
  sample_x[2] = -0.55158350755530561;
  sample_y[2] = 0.10316701511061116;

  for(int p = 0; p < 3; p++) {
    bool converged = false;
    for(int step = 0; step < 10; step++) {
      double surf = eval_at_pt(sample_x[p], sample_y[p], s_modal);
      double dsdx = eval_at_pt(sample_x[p], sample_y[p], dsdx_modal);
      double dsdy = eval_at_pt(sample_x[p], sample_y[p], dsdy_modal);

      double sqrnorm = dsdx * dsdx + dsdy * dsdy;
      if(sqrnorm > 0.0) {
        dsdx *= surf / sqrnorm;
        dsdy *= surf / sqrnorm;
      }

      sample_x[p] -= dsdx;
      sample_y[p] -= dsdy;

      // Converged
      if(dsdx * dsdx + dsdy * dsdy < *h * 1e-2) {
        converged = true;
        break;
      }
    }

    // Convert to x-y coords
    if(converged) {
      double r = sample_x[p];
      double s = sample_y[p];
      sample_x[p] = 0.5 * (-(r + s) * x[0] + (1.0 + r) * x[1] + (1.0 + s) * x[2]);
      sample_y[p] = 0.5 * (-(r + s) * y[0] + (1.0 + r) * y[1] + (1.0 + s) * y[2]);

      // Check within original element
      if(!is_point_in_cell(sample_x[p], sample_y[p], x, y)) {
        sample_x[p] = 0.0;
        sample_y[p] = 0.0;
      }
    } else {
      // TODO handle this properly
      sample_x[p] = 0.0;
      sample_y[p] = 0.0;
    }
  }
}
