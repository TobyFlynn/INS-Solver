inline void sigma_bflux(const int *bedgeNum, const double *sJ, const double *nx,
                        const double *ny, const double *s, const double *vis,
                        double *sigFx, double *sigFy) {
  int exInd = *bedgeNum * 6;

  double k = sqrt(*vis);

  for(int i = 0; i < 6; i++) {
    sigFx[exInd + i] += gaussW_g[i] * sJ[exInd + i] * nx[exInd + i] * k * s[exInd + i];
    sigFy[exInd + i] += gaussW_g[i] * sJ[exInd + i] * ny[exInd + i] * k * s[exInd + i];
  }
}
