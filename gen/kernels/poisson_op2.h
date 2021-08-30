inline void poisson_op2(const int *edgeNum, const bool *rev,
                        const double *mDL, const double *mDR,
                        const double *pDL, const double *pDR,
                        const double *gVPL, const double *gVPR,
                        const double *sJL, const double *sJR,
                        const double *hL, const double *hR,
                        const double *gFactorL, const double *gFactorR,
                        const double *factorL, const double *factorR,
                        double *op1L, double *op1R,
                        double *op2L, double *op2R) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  // Get right matrix for this edge
  const double *gVML, *gVMR;
  if(edgeL == 0) {
    gVML = gFInterp0_g;
  } else if(edgeL == 1) {
    gVML = gFInterp1_g;
  } else {
    gVML = gFInterp2_g;
  }

  if(edgeR == 0) {
    gVMR = gFInterp0_g;
  } else if(edgeR == 1) {
    gVMR = gFInterp1_g;
  } else {
    gVMR = gFInterp2_g;
  }

  // First edge term
  // gVM'*gw*rho^-1*gDnM
  // gVM'*gw*rho^-1*gDnP
  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      op2L[c_ind] = 0.0;
      op2R[c_ind] = 0.0;
      for(int k = 0; k < 6; k++) {
        // mD
        int b_ind = k * 10 + j;
        // Transpose of gVM
        int a_ind = k * 10 + i;
        // Rho and sJ ind
        int factors_indL = edgeL * 6 + k;
        int factors_indR = edgeR * 6 + k;
        int factors_indLR;
        int factors_indRR;
        if(reverse) {
          factors_indLR = edgeL * 6 + 6 - 1 - k;
          factors_indRR = edgeR * 6 + 6 - 1 - k;
        } else {
          factors_indLR = edgeL * 6 + k;
          factors_indRR = edgeR * 6 + k;
        }

        // op1L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
        //                * gFactorL[factors_indL] * mDL[b_ind];
        // op1R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
        //                * gFactorR[factors_indR] * mDR[b_ind];
        //
        // op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
        //                * gFactorR[factors_indRR] * pDL[b_ind];
        // op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
        //                * gFactorL[factors_indLR] * pDR[b_ind];

        op1L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * mDL[b_ind];
        op1R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * mDR[b_ind];

        op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * pDL[b_ind];
        op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * pDR[b_ind];
      }
    }
  }

  // Second edge term
  // rho^-1*gDnM'*gw*gVM
  // rho^-1*gDnM'*gw*gVP
  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      for(int k = 0; k < 6; k++) {
        // gVM and gVP
        int b_ind = k * 10 + j;
        // Transpose of mD
        int a_ind = k * 10 + i;
        // Rho and sJ ind
        int factors_indL = edgeL * 6 + k;
        int factors_indR = edgeR * 6 + k;
        int factors_indLR;
        int factors_indRR;
        if(reverse) {
          factors_indLR = edgeL * 6 + 6 - 1 - k;
          factors_indRR = edgeR * 6 + 6 - 1 - k;
        } else {
          factors_indLR = edgeL * 6 + k;
          factors_indRR = edgeR * 6 + k;
        }

        // op1L[c_ind] += -0.5 * gFactorL[factors_indL] * mDL[a_ind] * gaussW_g[k]
        //                * sJL[factors_indL] * gVML[b_ind];
        // op1R[c_ind] += -0.5 * gFactorR[factors_indR] * mDR[a_ind] * gaussW_g[k]
        //                * sJR[factors_indR] * gVMR[b_ind];
        //
        // op2L[c_ind] += 0.5 * gFactorL[factors_indL] * mDL[a_ind] * gaussW_g[k]
        //                * sJL[factors_indL] * gVPL[b_ind];
        // op2R[c_ind] += 0.5 * gFactorR[factors_indR] * mDR[a_ind] * gaussW_g[k]
        //                * sJR[factors_indR] * gVPR[b_ind];

        op1L[c_ind] += -0.5 * mDL[a_ind] * gaussW_g[k]
                       * sJL[factors_indL] * gVML[b_ind];
        op1R[c_ind] += -0.5 * mDR[a_ind] * gaussW_g[k]
                       * sJR[factors_indR] * gVMR[b_ind];

        op2L[c_ind] += 0.5 * mDL[a_ind] * gaussW_g[k]
                       * sJL[factors_indL] * gVPL[b_ind];
        op2R[c_ind] += 0.5 * mDR[a_ind] * gaussW_g[k]
                       * sJR[factors_indR] * gVPR[b_ind];
      }
    }
  }

  // Calculate penalty parameter
  double tau[6];
  double max = 0.0;
  for(int i = 0; i < 6; i++) {
    int indL = edgeL * 6 + i;
    int indR;
    if(reverse)
      indR = edgeR * 6 + 6 - 1 - i;
    else
      indR = edgeR * 6 + i;
    // tauL[i] = 0.5 * (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(*hL * gFactorL[indL], *hR * gFactorR[indR]);
    tau[i] = (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(*hL * gFactorL[indL], *hR * gFactorR[indR]);
    // if(max < tau[i]) {
    //   max = tau[i];
    // }
  }

  // for(int i = 0; i < 6; i++) {
  //   tau[i] = max;
  // }

  // Third edge term
  // gVM'*gw*tau*gVM
  // gVM'*gw*tau*gVP
  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      for(int k = 0; k < 6; k++) {
        // gVM and gVP
        int b_ind = k * 10 + j;
        // Transpose of gVM
        int a_ind = k * 10 + i;
        // sJ ind
        int factors_indL = edgeL * 6 + k;
        int factors_indR = edgeR * 6 + k;

        int tau_indL = k;
        int tau_indR;
        if(reverse) {
          tau_indR = 6 - k - 1;
        } else {
          tau_indR = k;
        }

        // op1L[c_ind] += gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
        //                * tau[tau_indL] * gVML[b_ind];
        // op1R[c_ind] += gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
        //                * tau[tau_indR] * gVMR[b_ind];
        //
        // op2L[c_ind] += -gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
        //                * tau[tau_indL] * gVPL[b_ind];
        // op2R[c_ind] += -gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
        //                * tau[tau_indR] * gVPR[b_ind];

        op1L[c_ind] += 0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * tau[tau_indL] * gVML[b_ind];
        op1R[c_ind] += 0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * tau[tau_indR] * gVMR[b_ind];

        op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * tau[tau_indL] * gVPL[b_ind];
        op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * tau[tau_indR] * gVPR[b_ind];
      }
    }
  }
}
