inline void init_gauss_grad(double *rx, double *sx, double *ry,  double *sy,
                            double *Dx0, double *Dy0, double *Dx1, double *Dy1,
                            double *Dx2, double *Dy2) {
  // J = -xs.*yr + xr.*ys
  double J[18];
  for(int i = 0; i < 18; i++) {
    J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
  }

  // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
  for(int i = 0; i < 18; i++) {
    double rx_n = sy[i] / J[i];
    double sx_n = -ry[i] / J[i];
    double ry_n = -sx[i] / J[i];
    double sy_n = rx[i] / J[i];
    rx[i] = rx_n;
    sx[i] = sx_n;
    ry[i] = ry_n;
    sy[i] = sy_n;
  }

  for(int j = 0; j < 10; j++) {
    for(int i = 0; i < 6; i++) {
      int ind = j * 6 + i;
      Dx0[ind] = rx[i] * gF0Dr_g[ind] + sx[i] * gF0Ds_g[ind];
      Dy0[ind] = ry[i] * gF0Dr_g[ind] + sy[i] * gF0Ds_g[ind];
    }
  }

  for(int j = 0; j < 10; j++) {
    for(int i = 0; i < 6; i++) {
      int ind = j * 6 + i;
      Dx1[ind] = rx[i + 6] * gF1Dr_g[ind] + sx[i + 6] * gF1Ds_g[ind];
      Dy1[ind] = ry[i + 6] * gF1Dr_g[ind] + sy[i + 6] * gF1Ds_g[ind];
    }
  }

  for(int j = 0; j < 10; j++) {
    for(int i = 0; i < 6; i++) {
      int ind = j * 6 + i;
      Dx2[ind] = rx[i + 2 * 6] * gF2Dr_g[ind] + sx[i + 2 * 6] * gF2Ds_g[ind];
      Dy2[ind] = ry[i + 2 * 6] * gF2Dr_g[ind] + sy[i + 2 * 6] * gF2Ds_g[ind];
    }
  }
}
