inline void poisson_matrix_2d_op2_partial(const int **p, const DG_FP *gF0Dr,
                    const DG_FP *gF0Ds, const DG_FP *gF1Dr, const DG_FP *gF1Ds,
                    const DG_FP *gF2Dr, const DG_FP *gF2Ds,
                    const DG_FP *gFInterp0, const DG_FP *gFInterp1,
                    const DG_FP *gFInterp2, const int *edgeNum,
                    const bool *reverse, const DG_FP **x, const DG_FP **y,
                    const DG_FP **sJ, const DG_FP **nx, const DG_FP **ny,
                    const DG_FP **h, DG_FP *op1L, DG_FP *op1R) {
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  // Get constants
  // Using same Gauss points so should be able to replace dg_gf_npL and
  // dg_gf_npR with DG_GF_NP
  const int dg_npL      = DG_CONSTANTS[(p[0][0] - 1) * DG_NUM_CONSTANTS];
  const int dg_npfL     = DG_CONSTANTS[(p[0][0] - 1) * DG_NUM_CONSTANTS + 1];
  const int dg_gf_npL   = DG_CONSTANTS[(p[0][0] - 1) * DG_NUM_CONSTANTS + 4];
  const DG_FP *gaussWL = &gaussW_g[(p[0][0] - 1) * DG_GF_NP];
  const int dg_npR      = DG_CONSTANTS[(p[1][0] - 1) * DG_NUM_CONSTANTS];
  const int dg_npfR     = DG_CONSTANTS[(p[1][0] - 1) * DG_NUM_CONSTANTS + 1];
  const int dg_gf_npR   = DG_CONSTANTS[(p[1][0] - 1) * DG_NUM_CONSTANTS + 4];
  const DG_FP *gaussWR = &gaussW_g[(p[1][0] - 1) * DG_GF_NP];

  const DG_FP *gDrL, *gDsL, *gDrR, *gDsR, *gVML, *gVMR;

  if(edgeL == 0) {
    gDrL = &gF0Dr[(p[0][0] - 1) * DG_GF_NP * DG_NP];
    gDsL = &gF0Ds[(p[0][0] - 1) * DG_GF_NP * DG_NP];
    gVML = &gFInterp0[(p[0][0] - 1) * DG_GF_NP * DG_NP];
  } else if(edgeL == 1) {
    gDrL = &gF1Dr[(p[0][0] - 1) * DG_GF_NP * DG_NP];
    gDsL = &gF1Ds[(p[0][0] - 1) * DG_GF_NP * DG_NP];
    gVML = &gFInterp1[(p[0][0] - 1) * DG_GF_NP * DG_NP];
  } else {
    gDrL = &gF2Dr[(p[0][0] - 1) * DG_GF_NP * DG_NP];
    gDsL = &gF2Ds[(p[0][0] - 1) * DG_GF_NP * DG_NP];
    gVML = &gFInterp2[(p[0][0] - 1) * DG_GF_NP * DG_NP];
  }

  if(edgeR == 0) {
    gDrR = &gF0Dr[(p[1][0] - 1) * DG_GF_NP * DG_NP];
    gDsR = &gF0Ds[(p[1][0] - 1) * DG_GF_NP * DG_NP];
    gVMR = &gFInterp0[(p[1][0] - 1) * DG_GF_NP * DG_NP];
  } else if(edgeR == 1) {
    gDrR = &gF1Dr[(p[1][0] - 1) * DG_GF_NP * DG_NP];
    gDsR = &gF1Ds[(p[1][0] - 1) * DG_GF_NP * DG_NP];
    gVMR = &gFInterp1[(p[1][0] - 1) * DG_GF_NP * DG_NP];
  } else {
    gDrR = &gF2Dr[(p[1][0] - 1) * DG_GF_NP * DG_NP];
    gDsR = &gF2Ds[(p[1][0] - 1) * DG_GF_NP * DG_NP];
    gVMR = &gFInterp2[(p[1][0] - 1) * DG_GF_NP * DG_NP];
  }

  DG_FP rxL[DG_GF_NP], sxL[DG_GF_NP], ryL[DG_GF_NP], syL[DG_GF_NP];
  DG_FP rxR[DG_GF_NP], sxR[DG_GF_NP], ryR[DG_GF_NP], syR[DG_GF_NP];

  // Left edge
  for(int m = 0; m < DG_GF_NP; m++) {
    rxL[m] = 0.0;
    sxL[m] = 0.0;
    ryL[m] = 0.0;
    syL[m] = 0.0;
    for(int n = 0; n < dg_npL; n++) {
      // int ind = m + n * DG_GF_NP;
      int ind = DG_MAT_IND(m, n, DG_GF_NP, dg_npL);
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
    for(int n = 0; n < dg_npR; n++) {
      // int ind = m + n * DG_GF_NP;
      int ind = DG_MAT_IND(m, n, DG_GF_NP, dg_npR);
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
  DG_FP mDL[DG_GF_NP * DG_NP], mDR[DG_GF_NP * DG_NP], pDL[DG_GF_NP * DG_NP], pDR[DG_GF_NP * DG_NP];
  for(int m = 0; m < DG_GF_NP; m++) {
    for(int n = 0; n < dg_npL; n++) {
      // int ind = m + n * DG_GF_NP;
      int ind = DG_MAT_IND(m, n, DG_GF_NP, dg_npL);
      int p_ind, p_norm_indR;
      if(*reverse) {
        // p_ind = DG_GF_NP - 1 - m + n * DG_GF_NP;
        p_ind = DG_MAT_IND(DG_GF_NP - 1 - m, n, DG_GF_NP, dg_npL);
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
    for(int n = 0; n < dg_npR; n++) {
      // int ind = m + n * DG_GF_NP;
      int ind = DG_MAT_IND(m, n, DG_GF_NP, dg_npR);
      int p_ind, p_norm_indL;
      if(*reverse) {
        // p_ind = DG_GF_NP - 1 - m + n * DG_GF_NP;
        p_ind = DG_MAT_IND(DG_GF_NP - 1 - m, n, DG_GF_NP, dg_npR);
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

  DG_FP max_hinv = fmax(h[0][edgeL * dg_npfL], h[1][edgeR * dg_npfR]);
  DG_FP tau = 0.5 * (DG_ORDER + 1) * (DG_ORDER + 2) * max_hinv;

  // Left edge
  for(int m = 0; m < dg_npL; m++) {
    for(int n = 0; n < dg_npL; n++) {
      // op col-major
      // int c_ind = m + n * dg_npL;
      int c_ind = DG_MAT_IND(m, n, dg_npL, dg_npL);
      for(int k = 0; k < DG_GF_NP; k++) {
        // Dx' and Dy'
        // int a_ind = m * DG_GF_NP + k;
        int a_ind = DG_MAT_IND(k, m, DG_GF_NP, dg_npL);
        // Dx and Dy
        // int b_ind = n * DG_GF_NP + k;
        int b_ind = DG_MAT_IND(k, n, DG_GF_NP, dg_npL);

        op1L[c_ind] += 0.5 * gaussWL[k] * sJ[0][exIndL + k] * tau * gVML[a_ind] * gVML[b_ind];
        op1L[c_ind] += -0.5 * gaussWL[k] * sJ[0][exIndL + k] * gVML[a_ind] * mDL[b_ind];
        op1L[c_ind] += -0.5 * gaussWL[k] * sJ[0][exIndL + k] * mDL[a_ind] * gVML[b_ind];
      }
    }
  }

  // Right edge
  for(int m = 0; m < dg_npR; m++) {
    for(int n = 0; n < dg_npR; n++) {
      // op col-major
      // int c_ind = m + n * dg_npR;
      int c_ind = DG_MAT_IND(m, n, dg_npR, dg_npR);
      for(int k = 0; k < DG_GF_NP; k++) {
        // Dx' and Dy'
        // int a_ind = m * DG_GF_NP + k;
        int a_ind = DG_MAT_IND(k, m, DG_GF_NP, dg_npR);
        // Dx and Dy
        // int b_ind = n * DG_GF_NP + k;
        int b_ind = DG_MAT_IND(k, n, DG_GF_NP, dg_npR);

        op1R[c_ind] += 0.5 * gaussWR[k] * sJ[1][exIndR + k] * tau * gVMR[a_ind] * gVMR[b_ind];
        op1R[c_ind] += -0.5 * gaussWR[k] * sJ[1][exIndR + k] * gVMR[a_ind] * mDR[b_ind];
        op1R[c_ind] += -0.5 * gaussWR[k] * sJ[1][exIndR + k] * mDR[a_ind] * gVMR[b_ind];
      }
    }
  }
}
