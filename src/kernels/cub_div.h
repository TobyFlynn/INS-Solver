inline void cub_div(const double *rx, const double *sx, const double *ry,
                    const double *sy, const double *J, double *temp0,
                    const double *temp1, const double *temp2, const double *temp3) {
  for(int i = 0; i < 46; i++) {
    double div = rx[i] * temp0[i] + sx[i] * temp1[i];
    div += ry[i] * temp2[i] + sy[i] * temp3[i];
    temp0[i] = cubW_g[i] * J[i] * div;
  }
}
