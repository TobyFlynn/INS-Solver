inline void cub_grad_w(double *temp0, const double *rx, const double *sx,
                     const double *ry, const double *sy, const double *J,
                     double *temp1, double *temp2, double *temp3) {
  for(int i = 0; i < 46; i++) {
    temp1[i] = cubW_g[i] * J[i] * sx[i] * temp0[i];
    temp2[i] = cubW_g[i] * J[i] * ry[i] * temp0[i];
    temp3[i] = cubW_g[i] * J[i] * sy[i] * temp0[i];
    temp0[i] = cubW_g[i] * J[i] * rx[i] * temp0[i];
  }
}
