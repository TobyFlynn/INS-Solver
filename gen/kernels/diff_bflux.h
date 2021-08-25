inline void diff_bflux(const int *bedgeNum, const double *sJ, const double *nx,
                       const double *ny, const double *s, const double *vis,
                       const double *sigX, const double *sigY, double *flux) {
  int exInd = *bedgeNum * 6;

  for(int i = 0; i < 6; i++) {
    flux[exInd + i] += gaussW_g[i] * sJ[exInd + i] * 0.5 * (*vis) * s[exInd + i];
    // flux[exInd + i] += gaussW_g[i] * sJ[exInd + i] * (0.5 * nx[exInd + i] * sigX[exInd + i] + 0.5 * ny[exInd + i] * sigY[exInd + i]);
    // flux[exInd + i] += fscale[exInd + i] * (nx[exInd + i] * sigX[fmask[i]] + ny[exInd + i] * sigY[fmask[i]]);
  }
}
