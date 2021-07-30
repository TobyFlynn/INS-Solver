inline void poisson_op5(const int *edgeType, const int *edgeNum,
                        const int *d0, const int *d1, const int *d2,
                        const double *mD0, const double *mD1,
                        const double *mD2, const double *sJ,
                        const double *h, const double *gFactor,
                        const double *factor, double *op) {
  // Get right matrices for this edge
  // (TODO switch matrices to be defined on edge set instead of cell set)
  const double *mD, *gVM;
  if(*edgeNum == 0) {
    mD  = mD0;
    gVM = gFInterp0_g;
  } else if(*edgeNum == 1) {
    mD  = mD1;
    gVM = gFInterp1_g;
  } else {
    mD  = mD2;
    gVM = gFInterp2_g;
  }

  for(int i = 0; i < 4 * 6; i++) {
    op[i] = 0.0;
  }

  if(*edgeType != *d0 && *edgeType != *d1 && *edgeType != *d2) {
    // First edge term
    // gVM'*gw*rho^-1*gDnM
    for(int i = 0; i < 4 * 6; i++) {
      int indT = (i % 4) * 6 + i / 4;
      int indSJ = *edgeNum * 4 + (i % 4);
      op[i] = gVM[indT] * gaussW_g[i % 4] * sJ[indSJ];
    }
  } else {
    // Calculate penalty parameter
    double tauA[4];
    for(int i = 0; i < 4; i++) {
      int ind = *edgeNum  * 4 + i;
      tauA[i] = 100 * 0.5 * 5 * 6 * (*h * gFactor[ind]);
      // tauA[i] = 100 * 0.5 * 5 * 6 * (*h);
    }
    // First edge term
    // gVM'*gw*rho^-1*gDnM
    for(int i = 0; i < 4 * 6; i++) {
      int indT = (i % 4) * 6 + i / 4;
      int indSJ = *edgeNum * 4 + (i % 4);
      int indFactor = (i / 4);

      op[i] = gVM[indT] * gaussW_g[i % 4] * sJ[indSJ] * tauA[i % 4]
              - factor[indFactor] * mD[indT] * gaussW_g[i % 4] * sJ[indSJ];
    }
  }
}
