inline void cub_grad(const double *rx, const double *sx,
                     const double *ry, const double *sy, const double *J,
                     double *temp0, double *temp1) {
  for(int i = 0; i < 46; i++) {
    double dru = temp0[i];
    double dsu = temp1[i];
    temp0[i] = cubW_g[i] * J[i] * (rx[i] * dru + sx[i] * dsu);
    temp1[i] = cubW_g[i] * J[i] * (ry[i] * dru + sy[i] * dsu);
  }
}
