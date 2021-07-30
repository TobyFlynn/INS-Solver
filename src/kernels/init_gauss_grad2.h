inline void init_gauss_grad2(const double *nx, const double *ny, const double *Dx0,
                             const double *Dy0, const double *Dx1, const double *Dy1,
                             const double *Dx2, const double *Dy2, double *d0,
                             double *d1, double *d2) {
  for(int m = 0; m < DG_GF_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind = m * DG_NP + n;
      d0[ind] = nx[m] * Dx0[ind] + ny[m] * Dy0[ind];
      d1[ind] = nx[m + DG_GF_NP] * Dx1[ind] + ny[m + DG_GF_NP] * Dy1[ind];
      d2[ind] = nx[m + 2 * DG_GF_NP] * Dx2[ind] + ny[m + 2 * DG_GF_NP] * Dy2[ind];
    }
  }
}
