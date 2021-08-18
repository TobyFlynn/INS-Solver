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

  for(int i = 0; i < DG_GF_NP * DG_NP; i++) {
    op[i] = 0.0;
  }

  if(*edgeType != *d0 && *edgeType != *d1 && *edgeType != *d2) {
    // First edge term
    // gVM'*gw*rho^-1*gDnM
    for(int i = 0; i < DG_GF_NP * DG_NP; i++) {
      int indT = (i % DG_GF_NP) * DG_NP + i / DG_GF_NP;
      int indSJ = *edgeNum * DG_GF_NP + (i % DG_GF_NP);
      op[i] = gVM[indT] * gaussW_g[i % DG_GF_NP] * sJ[indSJ];
    }
  } else {
    // Calculate penalty parameter
    double tauA[DG_GF_NP];
    double maxTau = 0.0;
    for(int i = 0; i < DG_GF_NP; i++) {
      int ind = *edgeNum  * DG_GF_NP + i;
      // tauA[i] = 0.5 * (DG_ORDER + 1) * (DG_ORDER + 2) * (*h * gFactor[ind]);
      tauA[i] = (DG_ORDER + 1) * (DG_ORDER + 2) * (*h * gFactor[ind]);
      // if(maxTau < tauA[i])
      //   maxTau = tauA[i];
    }

    // for(int i = 0; i < DG_GF_NP; i++) {
    //   tauA[i] = maxTau;
    // }

    // First edge term
    // gVM'*gw*rho^-1*gDnM
    for(int i = 0; i < DG_GF_NP * DG_NP; i++) {
      int indT = (i % DG_GF_NP) * DG_NP + i / DG_GF_NP;
      int indSJ = *edgeNum * DG_GF_NP + (i % DG_GF_NP);
      int indFactor = (i / DG_GF_NP);

      op[i] = gVM[indT] * gaussW_g[i % DG_GF_NP] * sJ[indSJ] * tauA[i % DG_GF_NP]
              - factor[indFactor] * mD[indT] * gaussW_g[i % DG_GF_NP] * sJ[indSJ];

      // op[i] = gVM[indT] * gaussW_g[i % DG_GF_NP] * sJ[indSJ] * tauA[i % DG_GF_NP]
      //         - factor[indSJ] * mD[indT] * gaussW_g[i % DG_GF_NP] * sJ[indSJ];
    }
  }
}
