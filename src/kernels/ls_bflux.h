inline void ls_bflux(const int *p, const int *bedgeNum, const double *sJ,
                     const double *nx, const double *ny, const double *s,
                     double *dsldx, double *dsrdx, double *dsldy,
                     double *dsrdy) {
  // Get constants
  const double *gaussW = &gaussW_g[(*p - 1) * DG_GF_NP];

  int exInd = *bedgeNum * DG_GF_NP;

  for(int i = 0; i < DG_GF_NP; i++) {
    dsldx[exInd + i] += gaussW[i] * sJ[exInd + i] * nx[exInd + i] * s[exInd + i];
    dsrdx[exInd + i] += gaussW[i] * sJ[exInd + i] * nx[exInd + i] * s[exInd + i];
    dsldy[exInd + i] += gaussW[i] * sJ[exInd + i] * ny[exInd + i] * s[exInd + i];
    dsrdy[exInd + i] += gaussW[i] * sJ[exInd + i] * ny[exInd + i] * s[exInd + i];
  }
}
