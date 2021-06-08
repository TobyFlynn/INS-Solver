inline void diff_bflux(const int *bedgeNum, const double *sJ, const double *nx,
                       const double *ny, const double *sigX, const double *sigY,
                       double *flux) {
  int exInd = 0;
  if(*bedgeNum == 1) exInd = 7;
  else if(*bedgeNum == 2) exInd = 2 * 7;

  for(int i = 0; i < 7; i++) {
    flux[exInd + i] += 0.5 * gaussW_g[i] * sJ[exInd + i] * (nx[exInd + i] * sigX[exInd + i] + ny[exInd + i] * sigY[exInd + i]);
    // flux[exInd + i] += fscale[exInd + i] * (nx[exInd + i] * sigX[fmask[i]] + ny[exInd + i] * sigY[fmask[i]]);
  }
}
