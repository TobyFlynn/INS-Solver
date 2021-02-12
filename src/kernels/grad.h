inline void grad(const double *div0, const double *div1, const double *rx,
                 const double *sx, const double *ry, const double *sy,
                 double *ux, double *uy) {
  for(int i = 0; i < 15; i++) {
    ux[i] = rx[i] * div0[i] + sx[i] * div1[i];
    uy[i] = ry[i] * div0[i] + sy[i] * div1[i];
  }
}
