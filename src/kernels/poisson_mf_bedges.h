inline void poisson_mf_bedges(const int *bedgeType, const int *bedgeNum,
                              const int *d0, const int *d1, const int *d2,
                              const double *sJ, const double *nx,
                              const double *ny, const double *tau,
                              const double *u, const double *dudx,
                              const double *dudy, double *fluxX, double *fluxY,
                              double *flux) {
  int exInd = 0;
  if(*bedgeNum == 1) exInd = 7;
  else if(*bedgeNum == 2) exInd = 2 * 7;

  if(*bedgeType == *d0 || *bedgeType == *d1 || *bedgeType == *d2) {
    for(int i = 0; i < 7; i++) {
      flux[exInd + i] += gaussW_g[i] * sJ[exInd + i] * (nx[exInd + i] * dudx[exInd + i] + ny[exInd + i] * dudy[exInd + i] - tau[*bedgeNum] * (u[exInd + i]));
    }
  } else {
    for(int i = 0; i < 7; i++) {
      fluxX[exInd + i] += nx[exInd + i] * gaussW_g[i] * sJ[exInd + i] * u[exInd + i];
      fluxY[exInd + i] += ny[exInd + i] * gaussW_g[i] * sJ[exInd + i] * u[exInd + i];
    }
  }
}
