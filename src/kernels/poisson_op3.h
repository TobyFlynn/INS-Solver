inline void poisson_op3(const int *edgeType, const int *edgeNum,
                        const int *d0, const int *d1, const int *d2,
                        const double *mD0, const double *mD1,
                        const double *mD2, const double *sJ,
                        const double *h, const double *gFactor,
                        const double *factor, double *op1) {
  if(*edgeType != *d0 && *edgeType != *d1 && *edgeType != *d2)
    return;
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

  // First edge term
  // gVM'*gw*rho^-1*gDnM
  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      for(int k = 0; k < 7; k++) {
        // mD
        int b_ind = k * 15 + j;
        // Transpose of gVM
        int ind = i * 7 + k;
        int a_ind = ((ind * 15) % (15 * 7)) + (ind / 7);
        // Rho and sJ ind
        int factors_ind = *edgeNum * 7 + k;

        op1[c_ind] += -0.5 * gVM[a_ind] * gaussW_g[k] * sJ[factors_ind]
                      * gFactor[factors_ind] * mD[b_ind];
      }
    }
  }

  // Second edge term
  // rho^-1*gDnM'*gw*gVM
  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      for(int k = 0; k < 7; k++) {
        // gVM and gVP
        int b_ind = k * 15 + j;
        // Transpose of mD
        int ind = i * 7 + k;
        int a_ind = ((ind * 15) % (15 * 7)) + (ind / 7);
        // Rho and sJ ind
        int factors_ind = *edgeNum * 7 + k;

        op1[c_ind] += -factor[i] * mD[a_ind] * gaussW_g[k]
                      * sJ[factors_ind] * gVM[b_ind];
      }
    }
  }


  // Calculate penalty parameter
  double tauA[7];
  for(int i = 0; i < 7; i++) {
    int ind = *edgeNum  * 7 + i;
    tauA[i] = 100 * 0.5 * 5 * 6 * (*h * gFactor[ind]);
    // tauA[i] = 100 * 0.5 * 5 * 6 * (*h);
  }


  // Third edge term
  // gVM'*gw*tau*gVM
  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      for(int k = 0; k < 7; k++) {
        // gVM and gVP
        int b_ind = k * 15 + j;
        // Transpose of gVM
        int ind = i * 7 + k;
        int a_ind = ((ind * 15) % (15 * 7)) + (ind / 7);
        // sJ ind
        int factors_ind = *edgeNum * 7 + k;

        op1[c_ind] += gVM[a_ind] * gaussW_g[k] * sJ[factors_ind]
                      * tauA[k] * gVM[b_ind];
      }
    }
  }
}
