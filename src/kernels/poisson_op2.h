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
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int c_ind = i * DG_NP + j;
      op2L[c_ind] = 0.0;
      op2R[c_ind] = 0.0;
      for(int k = 0; k < DG_GF_NP; k++) {
        // mD
        int b_ind = k * DG_NP + j;
        // Transpose of gVM
        int a_ind = k * DG_NP + i;
        // Rho and sJ ind
        int factors_indL = edgeL * DG_GF_NP + k;
        int factors_indR = edgeR * DG_GF_NP + k;
        int factors_indLR;
        int factors_indRR;
        if(reverse) {
          factors_indLR = edgeL * DG_GF_NP + DG_GF_NP - 1 - k;
          factors_indRR = edgeR * DG_GF_NP + DG_GF_NP - 1 - k;
        } else {
          factors_indLR = edgeL * DG_GF_NP + k;
          factors_indRR = edgeR * DG_GF_NP + k;
        }

        op1L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * gFactorL[factors_indL] * mDL[b_ind];
        op1R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * gFactorR[factors_indR] * mDR[b_ind];

        op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * gFactorR[factors_indRR] * pDL[b_ind];
        op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * gFactorL[factors_indLR] * pDR[b_ind];
      }
    }
  }

  // Second edge term
  // rho^-1*gDnM'*gw*gVM
  // rho^-1*gDnM'*gw*gVP
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int c_ind = i * DG_NP + j;
      for(int k = 0; k < DG_GF_NP; k++) {
        // gVM and gVP
        int b_ind = k * DG_NP + j;
        // Transpose of mD
        int a_ind = k * DG_NP + i;
        // Rho and sJ ind
        int factors_indL = edgeL * DG_GF_NP + k;
        int factors_indR = edgeR * DG_GF_NP + k;
        int factors_indLR;
        int factors_indRR;
        if(reverse) {
          factors_indLR = edgeL * DG_GF_NP + DG_GF_NP - 1 - k;
          factors_indRR = edgeR * DG_GF_NP + DG_GF_NP - 1 - k;
        } else {
          factors_indLR = edgeL * DG_GF_NP + k;
          factors_indRR = edgeR * DG_GF_NP + k;
        }

        // op1L[c_ind] += -factorL[i] * mDL[a_ind] * gaussW_g[k]
        //                * sJL[factors_indL] * gVML[b_ind];
        // op1R[c_ind] += -factorR[i] * mDR[a_ind] * gaussW_g[k]
        //                * sJR[factors_indR] * gVMR[b_ind];
        //
        // op2L[c_ind] += factorL[i] * mDL[a_ind] * gaussW_g[k]
        //                * sJL[factors_indL] * gVPL[b_ind];
        // op2R[c_ind] += factorR[i] * mDR[a_ind] * gaussW_g[k]
        //                * sJR[factors_indR] * gVPR[b_ind];

        // op1L[c_ind] += -gFactorL[factors_indL] * mDL[a_ind] * gaussW_g[k]
        //                * sJL[factors_indL] * gVML[b_ind];
        // op1R[c_ind] += -gFactorR[factors_indR] * mDR[a_ind] * gaussW_g[k]
        //                * sJR[factors_indR] * gVMR[b_ind];
        //
        // op2L[c_ind] += gFactorL[factors_indL] * mDL[a_ind] * gaussW_g[k]
        //                * sJL[factors_indL] * gVPL[b_ind];
        // op2R[c_ind] += gFactorR[factors_indR] * mDR[a_ind] * gaussW_g[k]
        //                * sJR[factors_indR] * gVPR[b_ind];

        op1L[c_ind] += -0.5 * gFactorL[factors_indL] * mDL[a_ind] * gaussW_g[k]
                       * sJL[factors_indL] * gVML[b_ind];
        op1R[c_ind] += -0.5 * gFactorR[factors_indR] * mDR[a_ind] * gaussW_g[k]
                       * sJR[factors_indR] * gVMR[b_ind];

        op2L[c_ind] += 0.5 * gFactorL[factors_indL] * mDL[a_ind] * gaussW_g[k]
                       * sJL[factors_indL] * gVPL[b_ind];
        op2R[c_ind] += 0.5 * gFactorR[factors_indR] * mDR[a_ind] * gaussW_g[k]
                       * sJR[factors_indR] * gVPR[b_ind];
      }
    }
  }

  // Calculate penalty parameter
  double tauL[DG_GF_NP];
  double tauR[DG_GF_NP];
  double maxL = 0.0;
  double maxR = 0.0;
  for(int i = 0; i < DG_GF_NP; i++) {
    int indL = edgeL * DG_GF_NP + i;
    int indR;
    if(reverse)
      indR = edgeR * DG_GF_NP + DG_GF_NP - 1 - i;
    else
      indR = edgeR * DG_GF_NP + i;
    // tauL[i] = 0.5 * (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(*hL * gFactorL[indL], *hR * gFactorR[indR]);
    tauL[i] = (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(*hL * gFactorL[indL], *hR * gFactorR[indR]);
    // if(maxL < tauL[i]) {
    //   maxL = tauL[i];
    // }
  }
  for(int i = 0; i < DG_GF_NP; i++) {
    int indL;
    int indR = edgeR * DG_GF_NP + i;
    if(reverse)
      indL = edgeL * DG_GF_NP + DG_GF_NP - 1 - i;
    else
      indL = edgeL * DG_GF_NP + i;
    // tauR[i] = 0.5 * (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(*hL * gFactorL[indL], *hR * gFactorR[indR]);
    tauR[i] = (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(*hL * gFactorL[indL], *hR * gFactorR[indR]);
    // if(maxR < tauR[i]) {
    //   maxR = tauR[i];
    // }
  }

  // for(int i = 0; i < DG_GF_NP; i++) {
  //   tauL[i] = maxL;
  //   tauR[i] = maxR;
  // }

  // Third edge term
  // gVM'*gw*tau*gVM
  // gVM'*gw*tau*gVP
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int c_ind = i * DG_NP + j;
      for(int k = 0; k < DG_GF_NP; k++) {
        // gVM and gVP
        int b_ind = k * DG_NP + j;
        // Transpose of gVM
        int a_ind = k * DG_NP + i;
        // sJ ind
        int factors_indL = edgeL * DG_GF_NP + k;
        int factors_indR = edgeR * DG_GF_NP + k;

        // op1L[c_ind] += gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
        //                * tauL[k] * gVML[b_ind];
        // op1R[c_ind] += gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
        //                * tauR[k] * gVMR[b_ind];
        //
        // op2L[c_ind] += -gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
        //                * tauL[k] * gVPL[b_ind];
        // op2R[c_ind] += -gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
        //                * tauR[k] * gVPR[b_ind];

        op1L[c_ind] += 0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * tauL[k] * gVML[b_ind];
        op1R[c_ind] += 0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * tauR[k] * gVMR[b_ind];

        op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * tauL[k] * gVPL[b_ind];
        op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * tauR[k] * gVPR[b_ind];
      }
    }
  }
}
