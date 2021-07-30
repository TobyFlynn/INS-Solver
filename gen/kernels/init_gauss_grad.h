inline void init_gauss_grad(double *rx, double *sx, double *ry,  double *sy,
                            double *Dx0, double *Dy0, double *Dx1, double *Dy1,
                            double *Dx2, double *Dy2) {
  // J = -xs.*yr + xr.*ys
  double J[9];
  for(int i = 0; i < 9; i++) {
    J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
  }

  // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
  for(int i = 0; i < 9; i++) {
    double rx_n = sy[i] / J[i];
    double sx_n = -ry[i] / J[i];
    double ry_n = -sx[i] / J[i];
    double sy_n = rx[i] / J[i];
    rx[i] = rx_n;
    sx[i] = sx_n;
    ry[i] = ry_n;
    sy[i] = sy_n;
  }

  for(int m = 0; m < 3; m++) {
    for(int n = 0; n < 3; n++) {
      int ind = m * 3 + n;
      Dx0[ind] = rx[m] * gF0Dr_g[ind] + sx[m] * gF0Ds_g[ind];
      Dy0[ind] = ry[m] * gF0Dr_g[ind] + sy[m] * gF0Ds_g[ind];
    }
  }

  for(int m = 0; m < 3; m++) {
    for(int n = 0; n < 3; n++) {
      int ind = m * 3 + n;
      Dx1[ind] = rx[m + 3] * gF1Dr_g[ind] + sx[m + 3] * gF1Ds_g[ind];
      Dy1[ind] = ry[m + 3] * gF1Dr_g[ind] + sy[m + 3] * gF1Ds_g[ind];
    }
  }

  for(int m = 0; m < 3; m++) {
    for(int n = 0; n < 3; n++) {
      int ind = m * 3 + n;
      Dx2[ind] = rx[m + 2 * 3] * gF2Dr_g[ind] + sx[m + 2 * 3] * gF2Ds_g[ind];
      Dy2[ind] = ry[m + 2 * 3] * gF2Dr_g[ind] + sy[m + 2 * 3] * gF2Ds_g[ind];
    }
  }
}
