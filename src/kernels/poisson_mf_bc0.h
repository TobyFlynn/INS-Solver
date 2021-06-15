inline void poisson_mf_bc0(const int *bedgeType, const int *bedgeNum,
                           const double *sJ, const double *nx, const double *ny,
                           const double *tau, const double *bc, double *fluxX,
                           double *fluxY, double *flux) {
  int exInd = 0;
  if(*bedgeNum == 1) exInd = 7;
  else if(*bedgeNum == 2) exInd = 2 * 7;

  for(int i = 0; i < 7; i++) {
    fluxX[exInd + i] += nx[exInd + i] * gaussW_g[i] * sJ[exInd + i] * bc[exInd + i];
    fluxY[exInd + i] += ny[exInd + i] * gaussW_g[i] * sJ[exInd + i] * bc[exInd + i];
    flux[exInd + i] += gaussW_g[i] * sJ[exInd + i] * tau[*bedgeNum] * bc[exInd + i];
  }
}
