inline void div(const double *div0, const double *div1, const double *div2,
                const double *div3, const double *rx, const double *sx,
                const double *ry, const double *sy, double *res) {
  for(int i = 0; i < 15; i++) {
    res[i] = rx[i] * div0[i] + sx[i] * div1[i] + ry[i] * div2[i] + sy[i] * div3[i];
  }
}
