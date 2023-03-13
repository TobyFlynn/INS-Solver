inline void fpmf_2d_apply_bc(const int *p, const DG_FP *gF0Dr,
                             const DG_FP *gF0Ds, const DG_FP *gF1Dr,
                             const DG_FP *gF1Ds, const DG_FP *gF2Dr,
                             const DG_FP *gF2Ds, const DG_FP *gFInterp0,
                             const DG_FP *gFInterp1, const DG_FP *gFInterp2,
                             const int *btype, const int *edgeNum,
                             const DG_FP *x, const DG_FP *y, const DG_FP *sJ,
                             const DG_FP *nx, const DG_FP *ny, const DG_FP *h,
                             const DG_FP *factor, const DG_FP *bc, DG_FP *rhs) {
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
  const int exInd = *edgeNum * dg_gf_np;

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
    DG_FP hinv = h[*edgeNum * dg_npf];
    for(int i = 0; i < DG_GF_NP; i++) {
      int ind = *edgeNum * DG_GF_NP + i;
      tau[i] = 0.5 * (DG_ORDER + 1) * (DG_ORDER + 2) * hinv * factor[ind];
    }

    // Apply BC matrix
    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_gf_np; j++) {
        int mat_ind = DG_MAT_IND(j, i, dg_gf_np, dg_np);
        rhs[i] += sJ[exInd + j] * gaussW[j] * gVM[mat_ind] * tau[j] * bc[j];
        rhs[i] -= sJ[exInd + j] * gaussW[j] * mD[mat_ind] * bc[j];
      }
    }
    /*
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
    */
  } else {
    // Neumann
    // Apply BC matrix
    // TODO double check
    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_gf_np; j++) {
        int mat_ind = DG_MAT_IND(j, i, dg_gf_np, dg_np);
        rhs[i] += sJ[exInd + j] * gaussW[j] * gVM[mat_ind] * bc[j];
      }
    }
    /*
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
    */
  }
}
