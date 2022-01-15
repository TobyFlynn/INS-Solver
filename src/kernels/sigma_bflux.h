inline void sigma_bflux(const int *p, const int *bedgeNum, const double *sJ,
                        const double *nx, const double *ny, const double *s,
                        const double *vis, double *sigFx, double *sigFy) {
  // Get constants
  const double *gaussW = &gaussW_g[(*p - 1) * DG_GF_NP];

  int exInd = *bedgeNum * DG_GF_NP;

  double k = sqrt(*vis);

  for(int i = 0; i < DG_GF_NP; i++) {
    sigFx[exInd + i] += gaussW[i] * sJ[exInd + i] * nx[exInd + i] * k * s[exInd + i];
    sigFy[exInd + i] += gaussW[i] * sJ[exInd + i] * ny[exInd + i] * k * s[exInd + i];
  }
}
