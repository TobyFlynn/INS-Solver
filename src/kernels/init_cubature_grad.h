inline void init_cubature_grad(double *rx, double *sx, double *ry,  double *sy,
                               double *Dx, double *Dy) {
  // J = -xs.*yr + xr.*ys
  double J[DG_CUB_NP];
  for(int i = 0; i < DG_CUB_NP; i++) {
    J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
  }

  // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
  for(int i = 0; i < DG_CUB_NP; i++) {
    double rx_n = sy[i] / J[i];
    double sx_n = -ry[i] / J[i];
    double ry_n = -sx[i] / J[i];
    double sy_n = rx[i] / J[i];
    rx[i] = rx_n;
    sx[i] = sx_n;
    ry[i] = ry_n;
    sy[i] = sy_n;
  }

  for(int j = 0; j < DG_NP; j++) {
    for(int i = 0; i < DG_CUB_NP; i++) {
      int ind = j * DG_CUB_NP + i;
      Dx[ind] = rx[i] * cubVDr_g[ind] + sx[i] * cubVDs_g[ind];
      Dy[ind] = ry[i] * cubVDr_g[ind] + sy[i] * cubVDs_g[ind];
    }
  }
}
