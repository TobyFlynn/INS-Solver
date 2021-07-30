inline void init_gauss_grad2(const double *nx, const double *ny, const double *Dx0,
                             const double *Dy0, const double *Dx1, const double *Dy1,
                             const double *Dx2, const double *Dy2, double *d0,
                             double *d1, double *d2) {
  for(int m = 0; m < 4; m++) {
    for(int n = 0; n < 6; n++) {
      int ind = m * 6 + n;
      d0[ind] = nx[m] * Dx0[ind] + ny[m] * Dy0[ind];
      d1[ind] = nx[m + 4] * Dx1[ind] + ny[m + 4] * Dy1[ind];
      d2[ind] = nx[m + 2 * 4] * Dx2[ind] + ny[m + 2 * 4] * Dy2[ind];
    }
  }
}
