inline void sigma_bflux(const int *bedgeNum, const double *sJ, const double *nx,
                        const double *ny, const double *s, double *sigFx,
                        double *sigFy) {
  int exInd = 0;
  if(*bedgeNum == 1) exInd = 7;
  else if(*bedgeNum == 2) exInd = 2 * 7;

  for(int i = 0; i < 7; i++) {
    sigFx[exInd + i] += gaussW_g[i] * sJ[exInd + i] * nx[exInd + i] * s[exInd + i];
    sigFy[exInd + i] += gaussW_g[i] * sJ[exInd + i] * ny[exInd + i] * s[exInd + i];
  }
}
