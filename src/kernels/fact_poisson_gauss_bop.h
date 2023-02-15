inline void fact_poisson_gauss_bop(const int *p, const DG_FP *gF0Dr,
                                   const DG_FP *gF0Ds, const DG_FP *gF1Dr,
                                   const DG_FP *gF1Ds, const DG_FP *gF2Dr,
                                   const DG_FP *gF2Ds, const DG_FP *gFInterp0,
                                   const DG_FP *gFInterp1, const DG_FP *gFInterp2,
                                   const int *btype, const int *edgeNum,
                                   const DG_FP *x, const DG_FP *y,
                                   const DG_FP *sJ, const DG_FP *nx,
                                   const DG_FP *ny, const DG_FP *h,
                                   const DG_FP *factor, DG_FP *op1,
                                   DG_FP *op_bc) {
  const DG_FP *gVM;
  if(*edgeNum == 0) {
    gVM = &gFInterp0[(*p - 1) * DG_GF_NP * DG_NP];
  } else if(*edgeNum == 1) {
    gVM = &gFInterp1[(*p - 1) * DG_GF_NP * DG_NP];
  } else {
    gVM = &gFInterp2[(*p - 1) * DG_GF_NP * DG_NP];
  }

  // Get constants
  const int dg_np      = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf     = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 1];
  const int dg_gf_np   = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 4];
  const DG_FP *gaussW = &gaussW_g[(*p - 1) * DG_GF_NP];

  // Dirichlet
  if(*btype == 0) {
    const DG_FP *gDr, *gDs;
    if(*edgeNum == 0) {
      gDr = &gF0Dr[(*p - 1) * DG_GF_NP * DG_NP];
      gDs = &gF0Ds[(*p - 1) * DG_GF_NP * DG_NP];
    } else if(*edgeNum == 1) {
      gDr = &gF1Dr[(*p - 1) * DG_GF_NP * DG_NP];
      gDs = &gF1Ds[(*p - 1) * DG_GF_NP * DG_NP];
    } else {
      gDr = &gF2Dr[(*p - 1) * DG_GF_NP * DG_NP];
      gDs = &gF2Ds[(*p - 1) * DG_GF_NP * DG_NP];
    }

    DG_FP rx[DG_GF_NP], sx[DG_GF_NP], ry[DG_GF_NP], sy[DG_GF_NP];

    for(int m = 0; m < dg_gf_np; m++) {
      rx[m] = 0.0;
      sx[m] = 0.0;
      ry[m] = 0.0;
      sy[m] = 0.0;
      for(int n = 0; n < dg_np; n++) {
        // int ind = m + n * dg_gf_np;
        int ind = DG_MAT_IND(m, n, dg_gf_np, dg_np);
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

    const int exInd = *edgeNum * dg_gf_np;
    DG_FP mD[DG_GF_NP * DG_NP];
    for(int m = 0; m < dg_gf_np; m++) {
      for(int n = 0; n < dg_np; n++) {
        // int ind = m + n * dg_gf_np;
        int ind = DG_MAT_IND(m, n, dg_gf_np, dg_np);

        DG_FP Dx = rx[m] * gDr[ind] + sx[m] * gDs[ind];
        DG_FP Dy = ry[m] * gDr[ind] + sy[m] * gDs[ind];
        mD[ind]   = factor[exInd + m] * (nx[exInd + m] * Dx + ny[exInd + m] * Dy);
      }
    }

    DG_FP tau[DG_GF_NP];
    for(int i = 0; i < DG_GF_NP; i++) {
      int ind = *edgeNum * DG_GF_NP + i;
      tau[i] = 0.5 * (*p + 1) * (*p + 2) * (*h) * factor[ind];
    }

    // Main matrix
    for(int m = 0; m < dg_np; m++) {
      for(int n = 0; n < dg_np; n++) {
        // op col-major
        // int c_ind = m + n * dg_np;
        int c_ind = DG_MAT_IND(m, n, dg_np, dg_np);
        // op row-major
        // int c_ind = m * dg_np + n;
        for(int k = 0; k < dg_gf_np; k++) {
          // Dx' and Dy'
          // int a_ind = m * dg_gf_np + k;
          int a_ind = DG_MAT_IND(k, m, dg_gf_np, dg_np);
          // Dx and Dy
          // int b_ind = n * dg_gf_np + k;
          int b_ind = DG_MAT_IND(k, n, dg_gf_np, dg_np);

          op1[c_ind] += gaussW[k] * sJ[exInd + k] * tau[k] * gVM[a_ind] * gVM[b_ind];
          op1[c_ind] += -gaussW[k] * sJ[exInd + k] * gVM[a_ind] * mD[b_ind];
          op1[c_ind] += -gaussW[k] * sJ[exInd + k] * mD[a_ind] * gVM[b_ind];
        }
      }
    }

    // Apply BC matrix
    for(int j = 0; j < dg_gf_np * dg_np; j++) {
      int indT_col = j;
      int col  = j % dg_gf_np;
      int row  = j / dg_gf_np;
      DG_FP val = gaussW[j % dg_gf_np] * sJ[*edgeNum * dg_gf_np + (j % dg_gf_np)] * tau[j % dg_gf_np];
      val *= gVM[indT_col];
      val -= mD[indT_col] * gaussW[j % dg_gf_np] * sJ[*edgeNum * dg_gf_np + (j % dg_gf_np)];
      int op_ind = DG_MAT_IND(row, col, dg_np, dg_gf_np);
      // op_bc[row + col * dg_np] = val;
      op_bc[op_ind] = val;
    }
  } else {
    // Neumann
    // Nothing for main matrix
    // Apply BC matrix
    for(int j = 0; j < dg_gf_np * dg_np; j++) {
      int indT_col = j;
      int col  = j % dg_gf_np;
      int row  = j / dg_gf_np;
      DG_FP val = gaussW[j % dg_gf_np] * sJ[*edgeNum * dg_gf_np + (j % dg_gf_np)];
      val *= gVM[indT_col];
      int op_ind = DG_MAT_IND(row, col, dg_np, dg_gf_np);
      // op_bc[row + col * dg_np] = val;
      op_bc[op_ind] = val;
    }
  }
}
