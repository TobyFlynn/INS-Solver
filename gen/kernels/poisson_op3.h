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
  for(int i = 0; i < 6; i++) {
    for(int j = 0; j < 6; j++) {
      int c_ind = i * 6 + j;
      for(int k = 0; k < 4; k++) {
        // mD
        int b_ind = k * 6 + j;
        // Transpose of gVM
        int ind = i * 4 + k;
        int a_ind = ((ind * 6) % (6 * 4)) + (ind / 4);
        // Rho and sJ ind
        int factors_ind = *edgeNum * 4 + k;

        op1[c_ind] += -0.5 * gVM[a_ind] * gaussW_g[k] * sJ[factors_ind]
                      * gFactor[factors_ind] * mD[b_ind];
      }
    }
  }

  // Second edge term
  // rho^-1*gDnM'*gw*gVM
  for(int i = 0; i < 6; i++) {
    for(int j = 0; j < 6; j++) {
      int c_ind = i * 6 + j;
      for(int k = 0; k < 4; k++) {
        // gVM and gVP
        int b_ind = k * 6 + j;
        // Transpose of mD
        int ind = i * 4 + k;
        int a_ind = ((ind * 6) % (6 * 4)) + (ind / 4);
        // Rho and sJ ind
        int factors_ind = *edgeNum * 4 + k;

        op1[c_ind] += -factor[i] * mD[a_ind] * gaussW_g[k]
                      * sJ[factors_ind] * gVM[b_ind];
      }
    }
  }


  // Calculate penalty parameter
  double tauA[4];
  for(int i = 0; i < 4; i++) {
    int ind = *edgeNum  * 4 + i;
    tauA[i] = 100 * 0.5 * 5 * 6 * (*h * gFactor[ind]);
    // tauA[i] = 100 * 0.5 * 5 * 6 * (*h);
  }


  // Third edge term
  // gVM'*gw*tau*gVM
  for(int i = 0; i < 6; i++) {
    for(int j = 0; j < 6; j++) {
      int c_ind = i * 6 + j;
      for(int k = 0; k < 4; k++) {
        // gVM and gVP
        int b_ind = k * 6 + j;
        // Transpose of gVM
        int ind = i * 4 + k;
        int a_ind = ((ind * 6) % (6 * 4)) + (ind / 4);
        // sJ ind
        int factors_ind = *edgeNum * 4 + k;

        op1[c_ind] += gVM[a_ind] * gaussW_g[k] * sJ[factors_ind]
                      * tauA[k] * gVM[b_ind];
      }
    }
  }
}
