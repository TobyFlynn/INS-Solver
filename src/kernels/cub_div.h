inline void cub_div(double *temp0, double *temp1, const double *rx, const double *sx,
                    const double *ry, const double *sy, const double *J,
                    double *temp2, double *temp3) {
  for(int i = 0; i < 46; i++) {
    double Vu = temp0[i];
    double Vv = temp1[i];
    temp0[i] = cubW[i] * J[i] * rx[i] * Vu;
    temp1[i] = cubW[i] * J[i] * sx[i] * Vu;
    temp2[i] = cubW[i] * J[i] * ry[i] * Vv;
    temp3[i] = cubW[i] * J[i] * sy[i] * Vv;
  }
}
