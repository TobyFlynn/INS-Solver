#include "utils.h"


inline void sample_interface(const double *r, const double *s,
                             const double *x, const double *y,
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
  if(!interface) {
    for(int i = 0; i < LS_SAMPLE_NP; i++) {
      sample_x[i] = NAN;
      sample_y[i] = NAN;
    }
    return;
  }

  // Initial positions of sample points (in r-s coords)
  for(int i = 0; i < LS_SAMPLE_NP; i++) {
    sample_x[i] = r[i] * 0.9;
    sample_y[i] = s[i] * 0.9;
  }

  for(int p = 0; p < LS_SAMPLE_NP; p++) {
    bool converged = false;
    for(int step = 0; step < 10; step++) {
      double surf = eval_at_pt(sample_x[p], sample_y[p], s_modal);
      double dsdx, dsdy;
      eval_grad_at_pt(sample_x[p], sample_y[p], s_modal, dsdx, dsdy);
      // double dsdx = eval_at_pt(sample_x[p], sample_y[p], dsdx_modal);
      // double dsdy = eval_at_pt(sample_x[p], sample_y[p], dsdy_modal);

      double sqrnorm = dsdx * dsdx + dsdy * dsdy;
      if(sqrnorm > 0.0) {
        dsdx *= surf / sqrnorm;
        dsdy *= surf / sqrnorm;
      }

      sample_x[p] -= dsdx;
      sample_y[p] -= dsdy;

      // Converged
      if(dsdx * dsdx + dsdy * dsdy < 1.5 * 1e-8) {
        converged = true;
        break;
      }
    }

    // Convert to x-y coords
    if(converged) {
      double r = sample_x[p];
      double s = sample_y[p];
      double new_x, new_y;
      rs_to_xy(r, s, new_x, new_y, x, y);
      sample_x[p] = new_x;
      sample_y[p] = new_y;

      // Check within original element
      if(!is_point_in_cell(sample_x[p], sample_y[p], x, y)) {
        sample_x[p] = NAN;
        sample_y[p] = NAN;
      }
    } else {
      // TODO handle this properly
      sample_x[p] = NAN;
      sample_y[p] = NAN;
    }
  }
}
