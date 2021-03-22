inline void init_gauss_grad(double *rx, double *sx, double *ry,  double *sy,
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

  for(int m = 0; m < 7; m++) {
    for(int n = 0; n < 15; n++) {
      int ind = m * 15 + n;
      Dx0[ind] = rx[m] * gF0Dr[ind] + sx[m] * gF0Ds[ind];
      Dy0[ind] = ry[m] * gF0Dr[ind] + sy[m] * gF0Ds[ind];
    }
  }

  for(int m = 0; m < 7; m++) {
    for(int n = 0; n < 15; n++) {
      int ind = m * 15 + n;
      Dx1[ind] = rx[m + 7] * gF1Dr[ind] + sx[m + 7] * gF1Ds[ind];
      Dy1[ind] = ry[m + 7] * gF1Dr[ind] + sy[m + 7] * gF1Ds[ind];
    }
  }

  for(int m = 0; m < 7; m++) {
    for(int n = 0; n < 15; n++) {
      int ind = m * 15 + n;
      Dx2[ind] = rx[m + 14] * gF2Dr[ind] + sx[m + 14] * gF2Ds[ind];
      Dy2[ind] = ry[m + 14] * gF2Dr[ind] + sy[m + 14] * gF2Ds[ind];
    }
  }
}
