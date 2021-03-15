inline void init_gauss_grad_neighbour(const int *reverse, double *rx,
                            double *sx, double *ry,  double *sy,
                            double *Dx0, double *Dy0, double *Dx1, double *Dy1,
                            double *Dx2, double *Dy2) {
  // J = -xs.*yr + xr.*ys
  double J[21];
  for(int i = 0; i < 21; i++) {
    J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
  }

  // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
  for(int i = 0; i < 21; i++) {
    double rx_n = sy[i] / J[i];
    double sx_n = -ry[i] / J[i];
    double ry_n = -sx[i] / J[i];
    double sy_n = rx[i] / J[i];
    rx[i] = rx_n;
    sx[i] = sx_n;
    ry[i] = ry_n;
    sy[i] = sy_n;
  }

  if(reverse[0]) {
    for(int m = 0; m < 7; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        Dx0[ind] = rx[m] * gF0DrR[ind] + sx[m] * gF0DsR[ind];
        Dy0[ind] = ry[m] * gF0DrR[ind] + sy[m] * gF0DsR[ind];
      }
    }
  } else {
    for(int m = 0; m < 7; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        Dx0[ind] = rx[m] * gF0Dr[ind] + sx[m] * gF0Ds[ind];
        Dy0[ind] = ry[m] * gF0Dr[ind] + sy[m] * gF0Ds[ind];
      }
    }
  }

  if(reverse[1]) {
    for(int m = 0; m < 7; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        Dx1[ind] = rx[m + 7] * gF1DrR[ind] + sx[m + 7] * gF1DsR[ind];
        Dy1[ind] = ry[m + 7] * gF1DrR[ind] + sy[m + 7] * gF1DsR[ind];
      }
    }
  } else {
    for(int m = 0; m < 7; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        Dx1[ind] = rx[m + 7] * gF1Dr[ind] + sx[m + 7] * gF1Ds[ind];
        Dy1[ind] = ry[m + 7] * gF1Dr[ind] + sy[m + 7] * gF1Ds[ind];
      }
    }
  }

  if(reverse[2]) {
    for(int m = 0; m < 7; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        Dx2[ind] = rx[m + 14] * gF2DrR[ind] + sx[m + 14] * gF2DsR[ind];
        Dy2[ind] = ry[m + 14] * gF2DrR[ind] + sy[m + 14] * gF2DsR[ind];
      }
    }
  } else {
    for(int m = 0; m < 7; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        Dx2[ind] = rx[m + 14] * gF2Dr[ind] + sx[m + 14] * gF2Ds[ind];
        Dy2[ind] = ry[m + 14] * gF2Dr[ind] + sy[m + 14] * gF2Ds[ind];
      }
    }
  }
}
