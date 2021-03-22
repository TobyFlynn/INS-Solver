inline void init_cubature(double *rx, double *sx, double *ry, double *sy,
                          double *J, double *temp) {
  // J = -xs.*yr + xr.*ys
  for(int i = 0; i < 46; i++) {
    J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
  }

  // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
  for(int i = 0; i < 46; i++) {
    double rx_n = sy[i] / J[i];
    double sx_n = -ry[i] / J[i];
    double ry_n = -sx[i] / J[i];
    double sy_n = rx[i] / J[i];
    rx[i] = rx_n;
    sx[i] = sx_n;
    ry[i] = ry_n;
    sy[i] = sy_n;
  }

  for(int m = 0; m < 46; m++) {
    for(int n = 0; n < 15; n++) {
      int ind = m * 15 + n;
      temp[ind] = J[m] * cubW[m] * cubV[ind];
    }
  }
}
