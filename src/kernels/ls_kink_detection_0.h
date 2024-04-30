inline void ls_kink_detection_0(const DG_FP *s, const DG_FP *dsdx, const DG_FP *dsdy, DG_FP **out) {
/*  
  const DG_FP d = s[0] / sqrt(dsdx[0] * dsdx[0] + dsdy[0] * dsdy[0]);
  const DG_FP n1_mag = sqrt(dsdx[1] * dsdx[1] + dsdy[1] * dsdy[1]);
  const DG_FP n1_x = dsdx[1] / n1_mag;
  const DG_FP n1_y = dsdy[1] / n1_mag;
  const DG_FP sqr_tol = 0.5 * 0.5;
  bool kink = false;
  for(int i = 2; i < DG_NP; i++) {
    const DG_FP n2_mag = sqrt(dsdx[i] * dsdx[i] + dsdy[i] * dsdy[i]);
    const DG_FP n2_x = dsdx[i] / n2_mag;
    const DG_FP n2_y = dsdy[i] / n2_mag;

    const DG_FP sqr_dist = (n1_x - n2_x) * (n1_x - n2_x) + (n1_y - n2_y) * (n1_y - n2_y);
    if(sqr_dist > sqr_tol) {
      kink = true;
      break;
    }
  }
*/

  const DG_FP sqr_tol = 0.5 * 0.5;
  // bool kink = false;
  for(int d_ind = 0; d_ind < DG_NP; d_ind++) {
    const DG_FP d = s[d_ind] / sqrt(dsdx[d_ind] * dsdx[d_ind] + dsdy[d_ind] * dsdy[d_ind]);
    for(int n1_ind = 0; n1_ind < DG_NP; n1_ind++) {
      if(n1_ind == d_ind) continue;
      const DG_FP n1_mag = sqrt(dsdx[n1_ind] * dsdx[n1_ind] + dsdy[n1_ind] * dsdy[n1_ind]);
      const DG_FP n1_x = dsdx[n1_ind] / n1_mag;
      const DG_FP n1_y = dsdy[n1_ind] / n1_mag;
      for(int n2_ind = 0; n2_ind < DG_NP; n2_ind++) {
        if(n2_ind == n1_ind || n2_ind == d_ind) continue;
        const DG_FP n2_mag = sqrt(dsdx[n2_ind] * dsdx[n2_ind] + dsdy[n2_ind] * dsdy[n2_ind]);
        const DG_FP n2_x = dsdx[n2_ind] / n2_mag;
        const DG_FP n2_y = dsdy[n2_ind] / n2_mag;

        const DG_FP sqr_dist = (n1_x - n2_x) * (n1_x - n2_x) + (n1_y - n2_y) * (n1_y - n2_y);
        if(sqr_dist > sqr_tol) {
          for(int i = 0; i < 3; i++) {
            out[i][0] += 1.0;
          }
          return;
        }
      }
    }
  }

/*
  for(int i = 0; i < DG_NP; i++) {
    DG_FP mag_sqr = dsdx[i] * dsdx[i] + dsdy[i] * dsdy[i];
    if(fabs(1.0 - mag_sqr) > 0.005) 
      out[i] = 1.0;
    else
      out[i] = 0.0;
  }
*/
}