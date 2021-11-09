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

  for(int j = 0; j < 10; j++) {
    for(int i = 0; i < 36; i++) {
      int ind = j * 36 + i;
      Dx[ind] = rx[i] * cubVDr_g[ind] + sx[i] * cubVDs_g[ind];
      Dy[ind] = ry[i] * cubVDr_g[ind] + sy[i] * cubVDs_g[ind];
    }
  }
}
