inline void pressure_solve_1(const int *edgeNum, const bool *rev,
                             const double **mD0, const double **mD1,
                             const double **mD2, const double **pD0,
                             const double **pD1, const double **pD2,
                             const double **gVP0, const double **gVP1,
                             const double **gVP2, const double **sJ,
                             const double **h, const double **tau, const double **gRho, const double **rho,
                             const double **u, double *rhsL, double *rhsR) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  // Get right matrices for this edge
  // (TODO switch matrices to be defined on edge set instead of cell set)
  const double *mDL, *mDR, *pDL, *pDR, *gVML, *gVMR, *gVPL, *gVPR;
  if(edgeL == 0) {
    mDL  = mD0[0];
    pDL  = pD0[0];
    gVML = gFInterp0_g;
    gVPL = gVP0[0];
  } else if(edgeL == 1) {
    mDL  = mD1[0];
    pDL  = pD1[0];
    gVML = gFInterp1_g;
    gVPL = gVP1[0];
  } else {
    mDL  = mD2[0];
    pDL  = pD2[0];
    gVML = gFInterp2_g;
    gVPL = gVP2[0];
  }

  if(edgeR == 0) {
    mDR  = mD0[1];
    pDR  = pD0[1];
    gVMR = gFInterp0_g;
    gVPR = gVP0[1];
  } else if(edgeR == 1) {
    mDR  = mD1[1];
    pDR  = pD1[1];
    gVMR = gFInterp1_g;
    gVPR = gVP1[1];
  } else {
    mDR  = mD2[1];
    pDR  = pD2[1];
    gVMR = gFInterp2_g;
    gVPR = gVP2[1];
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

        op1L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
                       * (1.0 / gRho[0][factors_indL]) * mDL[b_ind];
        op1R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
                       * (1.0 / gRho[1][factors_indR]) * mDR[b_ind];

        op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
                       * (1.0 / gRho[0][factors_indL]) * pDL[b_ind];
        op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
                       * (1.0 / gRho[1][factors_indR]) * pDR[b_ind];
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

        // op1L[c_ind] += -0.5 * (1.0 / rho[0][factors_indL]) * mDL[a_ind] * gaussW_g[k]
        //                * sJ[0][factors_indL] * gVML[b_ind];
        // op1R[c_ind] += -0.5 * (1.0 / rho[1][factors_indR]) * mDR[a_ind] * gaussW_g[k]
        //                * sJ[1][factors_indR] * gVMR[b_ind];
        //
        // op2L[c_ind] += 0.5 * (1.0 / rho[0][factors_indL]) * mDL[a_ind] * gaussW_g[k]
        //                * sJ[0][factors_indL] * gVPL[b_ind];
        // op2R[c_ind] += 0.5 * (1.0 / rho[1][factors_indR]) * mDR[a_ind] * gaussW_g[k]
        //                * sJ[1][factors_indR] * gVPR[b_ind];

        // op1L[c_ind] += -(1.0 / rho[0][factors_indL]) * mDL[a_ind] * gaussW_g[k]
        //                * sJ[0][factors_indL] * gVML[b_ind];
        // op1R[c_ind] += -(1.0 / rho[1][factors_indR]) * mDR[a_ind] * gaussW_g[k]
        //                * sJ[1][factors_indR] * gVMR[b_ind];
        //
        // op2L[c_ind] += (1.0 / rho[0][factors_indL]) * mDL[a_ind] * gaussW_g[k]
        //                * sJ[0][factors_indL] * gVPL[b_ind];
        // op2R[c_ind] += (1.0 / rho[1][factors_indR]) * mDR[a_ind] * gaussW_g[k]
        //                * sJ[1][factors_indR] * gVPR[b_ind];

        op1L[c_ind] += -(1.0 / rho[0][i]) * mDL[a_ind] * gaussW_g[k]
                       * sJ[0][factors_indL] * gVML[b_ind];
        op1R[c_ind] += -(1.0 / rho[1][i]) * mDR[a_ind] * gaussW_g[k]
                       * sJ[1][factors_indR] * gVMR[b_ind];

        op2L[c_ind] += (1.0 / rho[0][i]) * mDL[a_ind] * gaussW_g[k]
                       * sJ[0][factors_indL] * gVPL[b_ind];
        op2R[c_ind] += (1.0 / rho[1][i]) * mDR[a_ind] * gaussW_g[k]
                       * sJ[1][factors_indR] * gVPR[b_ind];
      }
    }
  }

  // Calculate penalty parameter
  double tauL[7];
  double tauR[7];
  for(int i = 0; i < 7; i++) {
    int indL = edgeL * 7 + i;
    int indR;
    if(reverse)
      indR = edgeR * 7 + 6 - i;
    else
      indR = edgeR * 7 + i;
    tauL[i] = 10 * 0.5 * 5 * 6 * fmax(*(h[0]) / gRho[0][indL], *(h[1]) / gRho[1][indR]);
  }
  for(int i = 0; i < 7; i++) {
    int indL;
    int indR = edgeR * 7 + i;
    if(reverse)
      indL = edgeL * 7 + 6 - i;
    else
      indL = edgeL * 7 + i;
    tauR[i] = 10 * 0.5 * 5 * 6 * fmax(*(h[0]) / gRho[0][indL], *(h[1]) / gRho[1][indR]);
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

        // op1L[c_ind] += 0.5 * gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
        //                * tau[0][edgeL] * gVML[b_ind];
        // op1R[c_ind] += 0.5 * gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
        //                * tau[1][edgeR] * gVMR[b_ind];
        //
        // op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
        //                * tau[0][edgeL] * gVPL[b_ind];
        // op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
        //                * tau[1][edgeR] * gVPR[b_ind];

        op1L[c_ind] += gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
                       * tauL[k] * gVML[b_ind];
        op1R[c_ind] += gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
                       * tauR[k] * gVMR[b_ind];

        op2L[c_ind] += -gVML[a_ind] * gaussW_g[k] * sJ[0][factors_indL]
                       * tauL[k] * gVPL[b_ind];
        op2R[c_ind] += -gVMR[a_ind] * gaussW_g[k] * sJ[1][factors_indR]
                       * tauR[k] * gVPR[b_ind];
      }
    }
  }


  // Multiply u by ops and add to rhs
  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int op_ind = i * 15 + j;
      rhsL[i] += op1L[op_ind] * u[0][j] + op2L[op_ind] * u[1][j];
      rhsR[i] += op1R[op_ind] * u[1][j] + op2R[op_ind] * u[0][j];
    }
  }
}
