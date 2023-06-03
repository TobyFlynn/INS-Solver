inline void poisson_coarse_matrix_over_int_2d_op2(const DG_FP *gF0Dr, const DG_FP *gF0Ds, 
                      const DG_FP *gF1Dr, const DG_FP *gF1Ds, const DG_FP *gF2Dr,
                      const DG_FP *gF2Ds, const DG_FP *gFInterp0, 
                      const DG_FP *gFInterp1, const DG_FP *gFInterp2, 
                      const int *edgeNum, const bool *reverse, const DG_FP **x, 
                      const DG_FP **y, const DG_FP **sJ, const DG_FP **nx,
                      const DG_FP **ny, const DG_FP **h, DG_FP *op1L, 
                      DG_FP *op1R, DG_FP *op2L, DG_FP *op2R) {
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  const DG_FP *gaussW = gaussW_g;

  const DG_FP *gDrL, *gDsL, *gDrR, *gDsR, *gVML, *gVMR;

  if(edgeL == 0) {
    gDrL = gF0Dr;
    gDsL = gF0Ds;
    gVML = gFInterp0;
  } else if(edgeL == 1) {
    gDrL = gF1Dr;
    gDsL = gF1Ds;
    gVML = gFInterp1;
  } else {
    gDrL = gF2Dr;
    gDsL = gF2Ds;
    gVML = gFInterp2;
  }

  if(edgeR == 0) {
    gDrR = gF0Dr;
    gDsR = gF0Ds;
    gVMR = gFInterp0;
  } else if(edgeR == 1) {
    gDrR = gF1Dr;
    gDsR = gF1Ds;
    gVMR = gFInterp1;
  } else {
    gDrR = gF2Dr;
    gDsR = gF2Ds;
    gVMR = gFInterp2;
  }

  DG_FP rxL[DG_GF_NP], sxL[DG_GF_NP], ryL[DG_GF_NP], syL[DG_GF_NP];
  DG_FP rxR[DG_GF_NP], sxR[DG_GF_NP], ryR[DG_GF_NP], syR[DG_GF_NP];

  // Left edge
  for(int m = 0; m < DG_GF_NP; m++) {
    rxL[m] = 0.0;
    sxL[m] = 0.0;
    ryL[m] = 0.0;
    syL[m] = 0.0;
    for(int n = 0; n < DG_NP_N1; n++) {
      // int ind = m + n * DG_GF_NP;
      int ind = DG_MAT_IND(m, n, DG_GF_NP, DG_NP_N1);
      rxL[m] += gDrL[ind] * x[0][n];
      sxL[m] += gDsL[ind] * x[0][n];
      ryL[m] += gDrL[ind] * y[0][n];
      syL[m] += gDsL[ind] * y[0][n];
    }
    DG_FP JL = -sxL[m] * ryL[m] + rxL[m] * syL[m];
    DG_FP rx_nL = syL[m] / JL;
    DG_FP sx_nL = -ryL[m] / JL;
    DG_FP ry_nL = -sxL[m] / JL;
    DG_FP sy_nL = rxL[m] / JL;
    rxL[m] = rx_nL;
    sxL[m] = sx_nL;
    ryL[m] = ry_nL;
    syL[m] = sy_nL;
  }

  // Right edge
  for(int m = 0; m < DG_GF_NP; m++) {
    rxR[m] = 0.0;
    sxR[m] = 0.0;
    ryR[m] = 0.0;
    syR[m] = 0.0;
    for(int n = 0; n < DG_NP_N1; n++) {
      // int ind = m + n * DG_GF_NP;
      int ind = DG_MAT_IND(m, n, DG_GF_NP, DG_NP_N1);
      rxR[m] += gDrR[ind] * x[1][n];
      sxR[m] += gDsR[ind] * x[1][n];
      ryR[m] += gDrR[ind] * y[1][n];
      syR[m] += gDsR[ind] * y[1][n];
    }
    DG_FP JR = -sxR[m] * ryR[m] + rxR[m] * syR[m];
    DG_FP rx_nR = syR[m] / JR;
    DG_FP sx_nR = -ryR[m] / JR;
    DG_FP ry_nR = -sxR[m] / JR;
    DG_FP sy_nR = rxR[m] / JR;
    rxR[m] = rx_nR;
    sxR[m] = sx_nR;
    ryR[m] = ry_nR;
    syR[m] = sy_nR;
  }

  // Left edge
  const int exIndL = edgeL * DG_GF_NP;
  const int exIndR = edgeR * DG_GF_NP;
  DG_FP mDL[DG_GF_NP * DG_NP_N1], mDR[DG_GF_NP * DG_NP_N1], pDL[DG_GF_NP * DG_NP_N1], pDR[DG_GF_NP * DG_NP_N1];
  for(int m = 0; m < DG_GF_NP; m++) {
    for(int n = 0; n < DG_NP_N1; n++) {
      // int ind = m + n * DG_GF_NP;
      int ind = DG_MAT_IND(m, n, DG_GF_NP, DG_NP_N1);
      int p_ind, p_norm_indR;
      if(*reverse) {
        // p_ind = DG_GF_NP - 1 - m + n * DG_GF_NP;
        p_ind = DG_MAT_IND(DG_GF_NP - 1 - m, n, DG_GF_NP, DG_NP_N1);
        p_norm_indR = exIndR + DG_GF_NP - 1 - m;
      } else {
        p_ind = ind;
        p_norm_indR = exIndR + m;
      }

      DG_FP DxL = rxL[m] * gDrL[ind] + sxL[m] * gDsL[ind];
      DG_FP DyL = ryL[m] * gDrL[ind] + syL[m] * gDsL[ind];
      mDL[ind]   = nx[0][exIndL + m] * DxL + ny[0][exIndL + m] * DyL;
      pDR[p_ind] = nx[1][p_norm_indR] * DxL + ny[1][p_norm_indR] * DyL;
    }
  }

  // Right edge
  for(int m = 0; m < DG_GF_NP; m++) {
    for(int n = 0; n < DG_NP_N1; n++) {
      // int ind = m + n * DG_GF_NP;
      int ind = DG_MAT_IND(m, n, DG_GF_NP, DG_NP_N1);
      int p_ind, p_norm_indL;
      if(*reverse) {
        // p_ind = DG_GF_NP - 1 - m + n * DG_GF_NP;
        p_ind = DG_MAT_IND(DG_GF_NP - 1 - m, n, DG_GF_NP, DG_NP_N1);
        p_norm_indL = exIndL + DG_GF_NP - 1 - m;
      } else {
        p_ind = ind;
        p_norm_indL = exIndL + m;
      }

      DG_FP DxR = rxR[m] * gDrR[ind] + sxR[m] * gDsR[ind];
      DG_FP DyR = ryR[m] * gDrR[ind] + syR[m] * gDsR[ind];
      mDR[ind]   = nx[1][exIndR + m] * DxR + ny[1][exIndR + m] * DyR;
      pDL[p_ind] = nx[0][p_norm_indL] * DxR + ny[0][p_norm_indL] * DyR;
    }
  }


  // Left edge
  DG_FP tauL[DG_GF_NP];
  DG_FP maxtau = 0.0;
  DG_FP max_hinv = fmax(h[0][edgeL * DG_NPF_N1], h[1][edgeR * DG_NPF_N1]);
  for(int i = 0; i < DG_GF_NP; i++) {
    int indL = edgeL * DG_GF_NP + i;
    int indR;
    if(reverse)
      indR = edgeR * DG_GF_NP + DG_GF_NP - 1 - i;
    else
      indR = edgeR * DG_GF_NP + i;

    tauL[i] = 0.5 * (DG_ORDER + 1) * (DG_ORDER + 2) * max_hinv;
    if(tauL[i] > maxtau) maxtau = tauL[i];
  }

  for(int i = 0; i < DG_GF_NP; i++) {
    tauL[i] = maxtau;
  }

  // Right edge
  DG_FP tauR[DG_GF_NP];
  maxtau = 0.0;
  for(int i = 0; i < DG_GF_NP; i++) {
    int indR = edgeR * DG_GF_NP + i;
    int indL;
    if(reverse)
      indL = edgeL * DG_GF_NP + DG_GF_NP - 1 - i;
    else
      indL = edgeL * DG_GF_NP + i;

    tauR[i] = 0.5 * (DG_ORDER + 1) * (DG_ORDER + 2) * max_hinv;
    if(tauR[i] > maxtau) maxtau = tauR[i];
  }

  for(int i = 0; i < DG_GF_NP; i++) {
    tauR[i] = maxtau;
  }

  // Left edge
  for(int m = 0; m < DG_NP_N1; m++) {
    for(int n = 0; n < DG_NP_N1; n++) {
      // op col-major
      // int c_ind = m + n * DG_NP_N1;
      int c_ind = DG_MAT_IND(m, n, DG_NP_N1, DG_NP_N1);
      for(int k = 0; k < DG_GF_NP; k++) {
        // Dx' and Dy'
        // int a_ind = m * DG_GF_NP + k;
        int a_ind = DG_MAT_IND(k, m, DG_GF_NP, DG_NP_N1);
        // Dx and Dy
        // int b_ind = n * DG_GF_NP + k;
        int b_ind = DG_MAT_IND(k, n, DG_GF_NP, DG_NP_N1);

        op1L[c_ind] += 0.5 * gaussW[k] * sJ[0][exIndL + k] * tauL[k] * gVML[a_ind] * gVML[b_ind];
        op1L[c_ind] += -0.5 * gaussW[k] * sJ[0][exIndL + k] * gVML[a_ind] * mDL[b_ind];
        op1L[c_ind] += -0.5 * gaussW[k] * sJ[0][exIndL + k] * mDL[a_ind] * gVML[b_ind];
      }
    }
  }

  for(int m = 0; m < DG_NP_N1; m++) {
    for(int n = 0; n < DG_NP_N1; n++) {
      // op col-major
      // int c_ind = m + n * DG_NP_N1;
      int c_ind = DG_MAT_IND(m, n, DG_NP_N1, DG_NP_N1);
      op2L[c_ind] = 0.0;
      for(int k = 0; k < DG_GF_NP; k++) {
        // Dx' and Dy'
        // int a_ind = m * DG_GF_NP + k;
        int a_ind = DG_MAT_IND(k, m, DG_GF_NP, DG_NP_N1);
        // Dx and Dy
        // int b_ind = n * DG_GF_NP + k;
        int b_ind = DG_MAT_IND(k, n, DG_GF_NP, DG_NP_N1);

        int b_indP;
        if(*reverse) {
          // b_indP = n * DG_GF_NP + DG_GF_NP - k - 1;
          b_indP = DG_MAT_IND(DG_GF_NP - k - 1, n, DG_GF_NP, DG_NP_N1);
        } else {
          // b_indP = n * DG_GF_NP + k;
          b_indP = DG_MAT_IND(k, n, DG_GF_NP, DG_NP_N1);
        }

        op2L[c_ind] += -0.5 * gaussW[k] * sJ[0][exIndL + k] * tauL[k] * gVML[a_ind] * gVMR[b_indP];
        op2L[c_ind] += -0.5 * gaussW[k] * sJ[0][exIndL + k] * gVML[a_ind] * pDL[b_ind];
        op2L[c_ind] += 0.5 * gaussW[k] * sJ[0][exIndL + k] * mDL[a_ind] * gVMR[b_indP];
      }
    }
  }

  // Right edge
  for(int m = 0; m < DG_NP_N1; m++) {
    for(int n = 0; n < DG_NP_N1; n++) {
      // op col-major
      // int c_ind = m + n * DG_NP_N1;
      int c_ind = DG_MAT_IND(m, n, DG_NP_N1, DG_NP_N1);
      for(int k = 0; k < DG_GF_NP; k++) {
        // Dx' and Dy'
        // int a_ind = m * DG_GF_NP + k;
        int a_ind = DG_MAT_IND(k, m, DG_GF_NP, DG_NP_N1);
        // Dx and Dy
        // int b_ind = n * DG_GF_NP + k;
        int b_ind = DG_MAT_IND(k, n, DG_GF_NP, DG_NP_N1);

        op1R[c_ind] += 0.5 * gaussW[k] * sJ[1][exIndR + k] * tauR[k] * gVMR[a_ind] * gVMR[b_ind];
        op1R[c_ind] += -0.5 * gaussW[k] * sJ[1][exIndR + k] * gVMR[a_ind] * mDR[b_ind];
        op1R[c_ind] += -0.5 * gaussW[k] * sJ[1][exIndR + k] * mDR[a_ind] * gVMR[b_ind];
      }
    }
  }

  for(int m = 0; m < DG_NP_N1; m++) {
    for(int n = 0; n < DG_NP_N1; n++) {
      // op col-major
      // int c_ind = m + n * DG_NP_N1;
      int c_ind = DG_MAT_IND(m, n, DG_NP_N1, DG_NP_N1);
      op2R[c_ind] = 0.0;
      for(int k = 0; k < DG_GF_NP; k++) {
        // Dx' and Dy'
        // int a_ind = m * DG_GF_NP + k;
        int a_ind = DG_MAT_IND(k, m, DG_GF_NP, DG_NP_N1);
        // Dx and Dy
        // int b_ind = n * DG_GF_NP + k;
        int b_ind = DG_MAT_IND(k, n, DG_GF_NP, DG_NP_N1);

        int b_indP;
        if(*reverse) {
          // b_indP = n * DG_GF_NP + DG_GF_NP - k - 1;
          b_indP = DG_MAT_IND(DG_GF_NP - k - 1, n, DG_GF_NP, DG_NP_N1);
        } else {
          // b_indP = n * DG_GF_NP + k;
          b_indP = DG_MAT_IND(k, n, DG_GF_NP, DG_NP_N1);
        }

        op2R[c_ind] += -0.5 * gaussW[k] * sJ[1][exIndR + k] * tauR[k] * gVMR[a_ind] * gVML[b_indP];
        op2R[c_ind] += -0.5 * gaussW[k] * sJ[1][exIndR + k] * gVMR[a_ind] * pDR[b_ind];
        op2R[c_ind] += 0.5 * gaussW[k] * sJ[1][exIndR + k] * mDR[a_ind] * gVML[b_indP];
      }
    }
  }
}
