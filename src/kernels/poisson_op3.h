inline void poisson_op3(const int *edgeType, const int *edgeNum,
                        const int *d0, const int *d1, const int *d2,
                        const double *mD, const double *sJ,
                        const double *h, const double *gFactor,
                        const double *factor, double *op1) {
  if(*edgeType != *d0 && *edgeType != *d1 && *edgeType != *d2)
    return;

  // Get right matrix for this edge
  const double *gVM;
  if(*edgeNum == 0) {
    gVM = gFInterp0_g;
  } else if(*edgeNum == 1) {
    gVM = gFInterp1_g;
  } else {
    gVM = gFInterp2_g;
  }

  // First edge term
  // gVM'*gw*rho^-1*gDnM
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int c_ind = i * DG_NP + j;
      for(int k = 0; k < DG_GF_NP; k++) {
        // mD
        int b_ind = k * DG_NP + j;
        // Transpose of gVM
        int ind = i * DG_GF_NP + k;
        int a_ind = ((ind * DG_NP) % (DG_NP * DG_GF_NP)) + (ind / DG_GF_NP);
        // Rho and sJ ind
        int factors_ind = *edgeNum * DG_GF_NP + k;

        // op1[c_ind] += -0.5 * gVM[a_ind] * gaussW_g[k] * sJ[factors_ind]
        //               * gFactor[factors_ind] * mD[b_ind];

        op1[c_ind] += -gVM[a_ind] * gaussW_g[k] * sJ[factors_ind]
                      * gFactor[factors_ind] * mD[b_ind];

        // op1[c_ind] += -gVM[a_ind] * gaussW_g[k] * sJ[factors_ind]
        //               * factor[i] * mD[b_ind];
      }
    }
  }

  // Second edge term
  // rho^-1*gDnM'*gw*gVM
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int c_ind = i * DG_NP + j;
      for(int k = 0; k < DG_GF_NP; k++) {
        // gVM and gVP
        int b_ind = k * DG_NP + j;
        // Transpose of mD
        int ind = i * DG_GF_NP + k;
        int a_ind = ((ind * DG_NP) % (DG_NP * DG_GF_NP)) + (ind / DG_GF_NP);
        // Rho and sJ ind
        int factors_ind = *edgeNum * DG_GF_NP + k;

        // op1[c_ind] += -factor[i] * mD[a_ind] * gaussW_g[k]
        //               * sJ[factors_ind] * gVM[b_ind];

        op1[c_ind] += -gFactor[factors_ind] * mD[a_ind] * gaussW_g[k]
                      * sJ[factors_ind] * gVM[b_ind];

        // op1[c_ind] += -0.5 * gFactor[factors_ind] * mD[a_ind] * gaussW_g[k]
        //               * sJ[factors_ind] * gVM[b_ind];

        // op1[c_ind] += -factor[i] * mD[a_ind] * gaussW_g[k]
        //               * sJ[factors_ind] * gVM[b_ind];
      }
    }
  }

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


  // Third edge term
  // gVM'*gw*tau*gVM
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int c_ind = i * DG_NP + j;
      for(int k = 0; k < DG_GF_NP; k++) {
        // gVM and gVP
        int b_ind = k * DG_NP + j;
        // Transpose of gVM
        int ind = i * DG_GF_NP + k;
        int a_ind = ((ind * DG_NP) % (DG_NP * DG_GF_NP)) + (ind / DG_GF_NP);
        // sJ ind
        int factors_ind = *edgeNum * DG_GF_NP + k;

        // op1[c_ind] += gVM[a_ind] * gaussW_g[k] * sJ[factors_ind]
        //               * tauA[k] * gVM[b_ind];

        op1[c_ind] += gVM[a_ind] * gaussW_g[k] * sJ[factors_ind]
                      * tauA[k] * gVM[b_ind];
      }
    }
  }
}
