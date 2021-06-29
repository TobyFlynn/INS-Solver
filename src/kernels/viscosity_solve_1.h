inline void viscosity_solve_1(const int *edgeNum, const bool *rev,
                              const double *mD0L, const double *mD0R,
                              const double *mD1L, const double *mD1R,
                              const double *mD2L, const double *mD2R,
                              const double *pD0L, const double *pD0R,
                              const double *pD1L, const double *pD1R,
                              const double *pD2L, const double *pD2R,
                              const double *gVP0L, const double *gVP0R,
                              const double *gVP1L, const double *gVP1R,
                              const double *gVP2L, const double *gVP2R,
                              const double *sJL, const double *sJR,
                              const double *hL, const double *hR,
                              const double *tauOL, const double *tauOR,
                              const double *gMuL, const double *gMuR,
                              const double *muL, const double *muR,
                              const double *rhoL, const double *rhoR,
                              const double *uL, const double *uR,
                              double *rhsL, double *rhsR) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  // Get right matrices for this edge
  // (TODO switch matrices to be defined on edge set instead of cell set)
  const double *mDL, *mDR, *pDL, *pDR, *gVML, *gVMR, *gVPL, *gVPR;
  if(edgeL == 0) {
    mDL  = mD0L;
    pDL  = pD0L;
    gVML = gFInterp0_g;
    gVPL = gVP0L;
  } else if(edgeL == 1) {
    mDL  = mD1L;
    pDL  = pD1L;
    gVML = gFInterp1_g;
    gVPL = gVP1L;
  } else {
    mDL  = mD2L;
    pDL  = pD2L;
    gVML = gFInterp2_g;
    gVPL = gVP2L;
  }

  if(edgeR == 0) {
    mDR  = mD0R;
    pDR  = pD0R;
    gVMR = gFInterp0_g;
    gVPR = gVP0R;
  } else if(edgeR == 1) {
    mDR  = mD1R;
    pDR  = pD1R;
    gVMR = gFInterp1_g;
    gVPR = gVP1R;
  } else {
    mDR  = mD2R;
    pDR  = pD2R;
    gVMR = gFInterp2_g;
    gVPR = gVP2R;
  }

  double op1L[15 * 15];
  double op1R[15 * 15];
  double op2L[15 * 15];
  double op2R[15 * 15];

  // First edge term
  // gVM'*gw*rho^-1*gDnM
  // gVM'*gw*rho^-1*gDnP
  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      op1L[c_ind] = 0.0;
      op1R[c_ind] = 0.0;
      op2L[c_ind] = 0.0;
      op2R[c_ind] = 0.0;
      for(int k = 0; k < 7; k++) {
        // mD
        int b_ind = k * 15 + j;
        // Transpose of gVM
        int ind = i * 7 + k;
        int a_ind = ((ind * 15) % (15 * 7)) + (ind / 7);
        // int a_ind = k * 15 + i;
        // Rho and sJ ind
        int factors_indL = edgeL * 7 + k;
        int factors_indR = edgeR * 7 + k;
        int factors_indLR;
        int factors_indRR;
        if(reverse) {
          factors_indLR = edgeL * 7 + 6 - k;
          factors_indRR = edgeR * 7 + 6 - k;
        } else {
          factors_indLR = edgeL * 7 + k;
          factors_indRR = edgeR * 7 + k;
        }
        /*
        op1L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
                       * gMu[0][factors_indL] * mDL[b_ind];
        op1R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
                       * gMu[1][factors_indR] * mDR[b_ind];

        op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
                       * gMu[1][factors_indRR] * pDL[b_ind];
        op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
                       * gMu[0][factors_indLR] * pDR[b_ind];
        */

        op1L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * gMuL[factors_indL] * mDL[b_ind];
        op1R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * gMuR[factors_indR] * mDR[b_ind];

        op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * gMuR[factors_indRR] * pDL[b_ind];
        op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * gMuL[factors_indLR] * pDR[b_ind];

        /*
        op1L[c_ind] += -0.5 * mu[0][i] * gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
                       * mDL[b_ind];
        op1R[c_ind] += -0.5 * mu[1][i] * gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
                       * mDR[b_ind];

        op2L[c_ind] += -0.5 * mu[0][i] * gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
                       * pDL[b_ind];
        op2R[c_ind] += -0.5 * mu[1][i] * gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
                       * pDR[b_ind];
        */
      }
    }
  }

  // Second edge term
  // rho^-1*gDnM'*gw*gVM
  // rho^-1*gDnM'*gw*gVP
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
        int factors_indL = edgeL * 7 + k;
        int factors_indR = edgeR * 7 + k;

        op1L[c_ind] += -muL[i] * mDL[a_ind] * gaussW_g[k]
                       * sJL[factors_indL] * gVML[b_ind];
        op1R[c_ind] += -muR[i] * mDR[a_ind] * gaussW_g[k]
                       * sJR[factors_indR] * gVMR[b_ind];

        op2L[c_ind] += muL[i] * mDL[a_ind] * gaussW_g[k]
                       * sJL[factors_indL] * gVPL[b_ind];
        op2R[c_ind] += muR[i] * mDR[a_ind] * gaussW_g[k]
                       * sJR[factors_indR] * gVPR[b_ind];

        /*
        op1L[c_ind] += -gMu[0][factors_indL] * mDL[a_ind] * gaussW_g[k]
                       * sJ[0][factors_indL] * gVML[b_ind];
        op1R[c_ind] += -gMu[1][factors_indR] * mDR[a_ind] * gaussW_g[k]
                       * sJ[1][factors_indR] * gVMR[b_ind];

        op2L[c_ind] += gMu[0][factors_indL] * mDL[a_ind] * gaussW_g[k]
                       * sJ[0][factors_indL] * gVPL[b_ind];
        op2R[c_ind] += gMu[1][factors_indR] * mDR[a_ind] * gaussW_g[k]
                       * sJ[1][factors_indR] * gVPR[b_ind];

        op1L[c_ind] += -mDLMu[a_ind] * gaussW_g[k]
                       * sJ[0][factors_indL] * gVML[b_ind];
        op1R[c_ind] += -mDRMu[a_ind] * gaussW_g[k]
                       * sJ[1][factors_indR] * gVMR[b_ind];

        op2L[c_ind] += mDLMu[a_ind] * gaussW_g[k]
                       * sJ[0][factors_indL] * gVPL[b_ind];
        op2R[c_ind] += mDRMu[a_ind] * gaussW_g[k]
                       * sJ[1][factors_indR] * gVPR[b_ind];

        op1L[c_ind] += -0.5 * mu[0][i] * mDL[a_ind] * gaussW_g[k]
                       * sJ[0][factors_indL] * gVML[b_ind];
        op1R[c_ind] += -0.5 * mu[1][i] * mDR[a_ind] * gaussW_g[k]
                       * sJ[1][factors_indR] * gVMR[b_ind];

        op2L[c_ind] += 0.5 * mu[0][i] * mDL[a_ind] * gaussW_g[k]
                       * sJ[0][factors_indL] * gVPL[b_ind];
        op2R[c_ind] += 0.5 * mu[1][i] * mDR[a_ind] * gaussW_g[k]
                       * sJ[1][factors_indR] * gVPR[b_ind];
        */
      }
    }
  }

  // Calculate penalty parameter
  double tauL[7];
  double tauR[7];
  double maxL = 0.0;
  double maxR = 0.0;
  for(int i = 0; i < 7; i++) {
    int indL = edgeL * 7 + i;
    int indR;
    if(reverse)
      indR = edgeR * 7 + 6 - i;
    else
      indR = edgeR * 7 + i;
    // tauL[i] = 10 * 0.5 * 5 * 6 * fmax(*hL / rhoL[indL], *hR / rhoR[indR]);
    tauL[i] = 10 * 0.5 * 5 * 6 * fmax(*hL * gMuL[indL], *hR * gMuR[indR]);
    if(maxL < tauL[i]) {
      maxL = tauL[i];
    }
  }
  for(int i = 0; i < 7; i++) {
    int indL;
    int indR = edgeR * 7 + i;
    if(reverse)
      indL = edgeL * 7 + 6 - i;
    else
      indL = edgeL * 7 + i;
    // tauR[i] = 10 * 0.5 * 5 * 6 * fmax(*hL / rhoL[indL], *hR / rhoR[indR]);
    tauR[i] = 10 * 0.5 * 5 * 6 * fmax(*hL * gMuL[indL], *hR * gMuR[indR]);
    if(maxR < tauR[i]) {
      maxR = tauR[i];
    }
  }

  for(int i = 0; i < 7; i++) {
    tauL[i] = maxL;
    tauR[i] = maxR;
  }

  // Third edge term
  // gVM'*gw*tau*gVM
  // gVM'*gw*tau*gVP
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
        int factors_indL = edgeL * 7 + k;
        int factors_indR = edgeR * 7 + k;

        op1L[c_ind] += gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * tauL[k] * gVML[b_ind];
        op1R[c_ind] += gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * tauR[k] * gVMR[b_ind];

        op2L[c_ind] += -gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * tauL[k] * gVPL[b_ind];
        op2R[c_ind] += -gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * tauR[k] * gVPR[b_ind];

        /*
        op1L[c_ind] += 0.5 * gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
                       * tau[0][edgeL] * gVML[b_ind];
        op1R[c_ind] += 0.5 * gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
                       * tau[1][edgeR] * gVMR[b_ind];

        op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
                       * tau[0][edgeL] * gVPL[b_ind];
        op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
                       * tau[1][edgeR] * gVPR[b_ind];
        */
      }
    }
  }


  // Multiply u by ops and add to rhs
  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int op_ind = i * 15 + j;
      rhsL[i] += op1L[op_ind] * uL[j] + op2L[op_ind] * uR[j];
      rhsR[i] += op1R[op_ind] * uR[j] + op2R[op_ind] * uL[j];
    }
  }
}
