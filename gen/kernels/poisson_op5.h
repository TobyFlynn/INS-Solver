inline void poisson_op5(const int *edgeType, const int *edgeNum,
                        const int *d0, const int *d1, const int *d2,
                        const double *mD, const double *sJ, const double *h,
                        const double *gFactor, const double *factor,
                        double *op) {
  // Get right matrix for this edge
  const double *gVM;
  if(*edgeNum == 0) {
    gVM = gFInterp0_g;
  } else if(*edgeNum == 1) {
    gVM = gFInterp1_g;
  } else {
    gVM = gFInterp2_g;
  }

  for(int i = 0; i < 6 * 10; i++) {
    op[i] = 0.0;
  }

  if(*edgeType != *d0 && *edgeType != *d1 && *edgeType != *d2) {
    // First edge term
    // gVM'*gw*rho^-1*gDnM
    for(int i = 0; i < 6 * 10; i++) {
      int indT = (i % 6) * 10 + i / 6;
      int indSJ = *edgeNum * 6 + (i % 6);
      op[i] = gVM[indT] * gaussW_g[i % 6] * sJ[indSJ];
    }
  } else {
    // Calculate penalty parameter
    double tauA[6];
    double maxTau = 0.0;
    for(int i = 0; i < 6; i++) {
      int ind = *edgeNum  * 6 + i;
      // tauA[i] = 0.5 * (DG_ORDER + 1) * (DG_ORDER + 2) * (*h * gFactor[ind]);
      tauA[i] = (DG_ORDER + 1) * (DG_ORDER + 2) * (*h * gFactor[ind]);
      // if(maxTau < tauA[i])
      //   maxTau = tauA[i];
    }

    // for(int i = 0; i < 6; i++) {
    //   tauA[i] = maxTau;
    // }

    // First edge term
    // gVM'*gw*rho^-1*gDnM
    for(int i = 0; i < 6 * 10; i++) {
      int indT = (i % 6) * 10 + i / 6;
      int indSJ = *edgeNum * 6 + (i % 6);
      int indFactor = (i / 6);

      op[i] = gVM[indT] * gaussW_g[i % 6] * sJ[indSJ] * tauA[i % 6]
              - factor[indFactor] * mD[indT] * gaussW_g[i % 6] * sJ[indSJ];

      // op[i] = gVM[indT] * gaussW_g[i % 6] * sJ[indSJ] * tauA[i % 6]
      //         - factor[indSJ] * mD[indT] * gaussW_g[i % 6] * sJ[indSJ];
    }
  }
}
