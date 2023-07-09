inline void factor_poisson_coarse_matrix_2d_bop_over_int(const DG_FP *gF0Dr, const DG_FP *gF0Ds, 
                            const DG_FP *gF1Dr, const DG_FP *gF1Ds, const DG_FP *gF2Dr,
                            const DG_FP *gF2Ds, const DG_FP *gFInterp0,
                            const DG_FP *gFInterp1, const DG_FP *gFInterp2,
                            const int *btype, const int *edgeNum, const DG_FP *x,
                            const DG_FP *y, const DG_FP *sJ, const DG_FP *nx,
                            const DG_FP *ny, const DG_FP *h, const DG_FP *factor,
                            DG_FP *op1, DG_FP *op_bc) {
  const DG_FP *gVM;
  if(*edgeNum == 0) {
    gVM = gFInterp0;
  } else if(*edgeNum == 1) {
    gVM = gFInterp1;
  } else {
    gVM = gFInterp2;
  }

  // Get constants
  const DG_FP *gaussW = gaussW_g;

  // Dirichlet
  if(*btype == 0) {
    const DG_FP *gDr, *gDs;
    if(*edgeNum == 0) {
      gDr = gF0Dr;
      gDs = gF0Ds;
    } else if(*edgeNum == 1) {
      gDr = gF1Dr;
      gDs = gF1Ds;
    } else {
      gDr = gF2Dr;
      gDs = gF2Ds;
    }

    DG_FP rx[DG_GF_NP], sx[DG_GF_NP], ry[DG_GF_NP], sy[DG_GF_NP];

    for(int m = 0; m < DG_GF_NP; m++) {
      rx[m] = 0.0;
      sx[m] = 0.0;
      ry[m] = 0.0;
      sy[m] = 0.0;
      for(int n = 0; n < DG_NP_N1; n++) {
        // int ind = m + n * DG_GF_NP;
        int ind = DG_MAT_IND(m, n, DG_GF_NP, DG_NP_N1);
        rx[m] += gDr[ind] * x[n];
        sx[m] += gDs[ind] * x[n];
        ry[m] += gDr[ind] * y[n];
        sy[m] += gDs[ind] * y[n];
      }
      DG_FP J = -sx[m] * ry[m] + rx[m] * sy[m];
      DG_FP rx_n = sy[m] / J;
      DG_FP sx_n = -ry[m] / J;
      DG_FP ry_n = -sx[m] / J;
      DG_FP sy_n = rx[m] / J;
      rx[m] = rx_n;
      sx[m] = sx_n;
      ry[m] = ry_n;
      sy[m] = sy_n;
    }

    const int exInd = *edgeNum * DG_GF_NP;
    DG_FP mD[DG_GF_NP * DG_NP_N1];
    for(int m = 0; m < DG_GF_NP; m++) {
      for(int n = 0; n < DG_NP_N1; n++) {
        // int ind = m + n * DG_GF_NP;
        int ind = DG_MAT_IND(m, n, DG_GF_NP, DG_NP_N1);

        DG_FP Dx = rx[m] * gDr[ind] + sx[m] * gDs[ind];
        DG_FP Dy = ry[m] * gDr[ind] + sy[m] * gDs[ind];
        mD[ind]   = factor[exInd + m] * (nx[exInd + m] * Dx + ny[exInd + m] * Dy);
      }
    }

    DG_FP tau[DG_GF_NP];
    DG_FP hinv = h[*edgeNum * DG_NPF_N1];
    for(int i = 0; i < DG_GF_NP; i++) {
      int ind = *edgeNum * DG_GF_NP + i;
      tau[i] = 0.5 * (DG_ORDER + 1) * (DG_ORDER + 2) * hinv * factor[ind];
    }

    // Main matrix
    for(int m = 0; m < DG_NP_N1; m++) {
      for(int n = 0; n < DG_NP_N1; n++) {
        // op col-major
        // int c_ind = m + n * DG_NP_N1;
        int c_ind = DG_MAT_IND(m, n, DG_NP_N1, DG_NP_N1);
        // op row-major
        // int c_ind = m * DG_NP_N1 + n;
        for(int k = 0; k < DG_GF_NP; k++) {
          // Dx' and Dy'
          // int a_ind = m * DG_GF_NP + k;
          int a_ind = DG_MAT_IND(k, m, DG_GF_NP, DG_NP_N1);
          // Dx and Dy
          // int b_ind = n * DG_GF_NP + k;
          int b_ind = DG_MAT_IND(k, n, DG_GF_NP, DG_NP_N1);

          op1[c_ind] += gaussW[k] * sJ[exInd + k] * tau[k] * gVM[a_ind] * gVM[b_ind];
          op1[c_ind] += -gaussW[k] * sJ[exInd + k] * gVM[a_ind] * mD[b_ind];
          op1[c_ind] += -gaussW[k] * sJ[exInd + k] * mD[a_ind] * gVM[b_ind];
        }
      }
    }

    // Apply BC matrix
    for(int j = 0; j < DG_GF_NP * DG_NP_N1; j++) {
      int indT_col = j;
      int col  = j % DG_GF_NP;
      int row  = j / DG_GF_NP;
      DG_FP val = gaussW[j % DG_GF_NP] * sJ[*edgeNum * DG_GF_NP + (j % DG_GF_NP)] * tau[j % DG_GF_NP];
      val *= gVM[indT_col];
      val -= mD[indT_col] * gaussW[j % DG_GF_NP] * sJ[*edgeNum * DG_GF_NP + (j % DG_GF_NP)];
      int op_ind = DG_MAT_IND(row, col, DG_NP_N1, DG_GF_NP);
      // op_bc[row + col * DG_NP_N1] = val;
      op_bc[op_ind] = val;
    }
  } else {
    // Neumann
    // Nothing for main matrix
    // Apply BC matrix
    for(int j = 0; j < DG_GF_NP * DG_NP_N1; j++) {
      int indT_col = j;
      int col  = j % DG_GF_NP;
      int row  = j / DG_GF_NP;
      DG_FP val = gaussW[j % DG_GF_NP] * sJ[*edgeNum * DG_GF_NP + (j % DG_GF_NP)];
      val *= gVM[indT_col];
      int op_ind = DG_MAT_IND(row, col, DG_NP_N1, DG_GF_NP);
      // op_bc[row + col * DG_NP_N1] = val;
      op_bc[op_ind] = val;
    }
  }
}
