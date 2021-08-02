inline void init_cubature_grad(double *rx, double *sx, double *ry,  double *sy,
                               double *Dx, double *Dy) {
  // J = -xs.*yr + xr.*ys
  double J[36];
  for(int i = 0; i < 36; i++) {
    J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
  }

  // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
  for(int i = 0; i < 36; i++) {
    double rx_n = sy[i] / J[i];
    double sx_n = -ry[i] / J[i];
    double ry_n = -sx[i] / J[i];
    double sy_n = rx[i] / J[i];
    rx[i] = rx_n;
    sx[i] = sx_n;
    ry[i] = ry_n;
    sy[i] = sy_n;
  }

  for(int m = 0; m < 36; m++) {
    for(int n = 0; n < 10; n++) {
      int ind = m * 10 + n;
      Dx[ind] = rx[m] * cubVDr_g[ind] + sx[m] * cubVDs_g[ind];
      Dy[ind] = ry[m] * cubVDr_g[ind] + sy[m] * cubVDs_g[ind];
    }
  }
}
