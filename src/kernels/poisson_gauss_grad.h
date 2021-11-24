inline void poisson_gauss_grad(const int **p, const double *gF0Dr,
                               const double *gF0Ds, const double *gF1Dr,
                               const double *gF1Ds, const double *gF2Dr,
                               const double *gF2Ds, const double *gFInterp0,
                               const double *gFInterp1, const double *gFInterp2,
                               const int *edgeNum, const bool *reverse,
                               const double **x, const double **y,
                               const double **sJ, const double **nx,
                               const double **ny, const double **fscale,
                               double *op1L, double *op1R, double *op2L,
                               double *op2R) {
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  // Get constants
  const int dg_np      = DG_CONSTANTS[(p[0][0] - 1) * 5];
  const int dg_npf     = DG_CONSTANTS[(p[0][0] - 1) * 5 + 1];
  const int dg_gf_np   = DG_CONSTANTS[(p[0][0] - 1) * 5 + 4];
  const double *gaussW = &gaussW_g[(p[0][0] - 1) * DG_GF_NP];

  const double *gDrL, *gDsL, *gDrR, *gDsR, *gVML, *gVMR;

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
    gDrR = &gF0Dr[(p[0][0] - 1) * DG_GF_NP * DG_NP];
    gDsR = &gF0Ds[(p[0][0] - 1) * DG_GF_NP * DG_NP];
    gVMR = &gFInterp0[(p[0][0] - 1) * DG_GF_NP * DG_NP];
  } else if(edgeR == 1) {
    gDrR = &gF1Dr[(p[0][0] - 1) * DG_GF_NP * DG_NP];
    gDsR = &gF1Ds[(p[0][0] - 1) * DG_GF_NP * DG_NP];
    gVMR = &gFInterp1[(p[0][0] - 1) * DG_GF_NP * DG_NP];
  } else {
    gDrR = &gF2Dr[(p[0][0] - 1) * DG_GF_NP * DG_NP];
    gDsR = &gF2Ds[(p[0][0] - 1) * DG_GF_NP * DG_NP];
    gVMR = &gFInterp2[(p[0][0] - 1) * DG_GF_NP * DG_NP];
  }

  double rxL[DG_GF_NP], sxL[DG_GF_NP], ryL[DG_GF_NP], syL[DG_GF_NP];
  double rxR[DG_GF_NP], sxR[DG_GF_NP], ryR[DG_GF_NP], syR[DG_GF_NP];

  for(int m = 0; m < dg_gf_np; m++) {
    rxL[m] = 0.0;
    sxL[m] = 0.0;
    ryL[m] = 0.0;
    syL[m] = 0.0;
    rxR[m] = 0.0;
    sxR[m] = 0.0;
    ryR[m] = 0.0;
    syR[m] = 0.0;
    for(int n = 0; n < dg_np; n++) {
      int ind = m + n * dg_gf_np;
      rxL[m] += gDrL[ind] * x[0][n];
      sxL[m] += gDsL[ind] * x[0][n];
      ryL[m] += gDrL[ind] * y[0][n];
      syL[m] += gDsL[ind] * y[0][n];
      rxR[m] += gDrR[ind] * x[1][n];
      sxR[m] += gDsR[ind] * x[1][n];
      ryR[m] += gDrR[ind] * y[1][n];
      syR[m] += gDsR[ind] * y[1][n];
    }
    double JL = -sxL[m] * ryL[m] + rxL[m] * syL[m];
    double JR = -sxR[m] * ryR[m] + rxR[m] * syR[m];
    double rx_nL = syL[m] / JL;
    double sx_nL = -ryL[m] / JL;
    double ry_nL = -sxL[m] / JL;
    double sy_nL = rxL[m] / JL;
    double rx_nR = syR[m] / JR;
    double sx_nR = -ryR[m] / JR;
    double ry_nR = -sxR[m] / JR;
    double sy_nR = rxR[m] / JR;
    rxL[m] = rx_nL;
    sxL[m] = sx_nL;
    ryL[m] = ry_nL;
    syL[m] = sy_nL;
    rxR[m] = rx_nR;
    sxR[m] = sx_nR;
    ryR[m] = ry_nR;
    syR[m] = sy_nR;
  }

  const int exIndL = edgeL * dg_gf_np;
  const int exIndR = edgeR * dg_gf_np;
  double mDL[DG_GF_NP * DG_NP], mDR[DG_GF_NP * DG_NP], pDL[DG_GF_NP * DG_NP], pDR[DG_GF_NP * DG_NP];
  for(int m = 0; m < dg_gf_np; m++) {
    for(int n = 0; n < dg_np; n++) {
      int ind = m + n * dg_gf_np;
      int p_ind, p_norm_indL, p_norm_indR;
      if(*reverse) {
        p_ind = dg_gf_np - 1 - m + n * dg_gf_np;
        p_norm_indL = exIndL + dg_gf_np - 1 - m;
        p_norm_indR = exIndR + dg_gf_np - 1 - m;
      } else {
        p_ind = ind;
        p_norm_indL = exIndL + m;
        p_norm_indR = exIndR + m;
      }

      double DxL = rxL[m] * gDrL[ind] + sxL[m] * gDsL[ind];
      double DyL = ryL[m] * gDrL[ind] + syL[m] * gDsL[ind];
      double DxR = rxR[m] * gDrR[ind] + sxR[m] * gDsR[ind];
      double DyR = ryR[m] * gDrR[ind] + syR[m] * gDsR[ind];
      mDL[ind]   = nx[0][exIndL + m] * DxL + ny[0][exIndL + m] * DyL;
      mDR[ind]   = nx[1][exIndR + m] * DxR + ny[1][exIndR + m] * DyR;
      pDL[p_ind] = nx[0][p_norm_indL] * DxR + ny[0][p_norm_indL] * DyR;
      pDR[p_ind] = nx[1][p_norm_indR] * DxL + ny[1][p_norm_indR] * DyL;
    }
  }

  double tau;
  if(fscale[0][edgeL * dg_npf] > fscale[1][edgeR * dg_npf]) {
    tau = 20 * 25 * fscale[0][edgeL * dg_npf];
  } else {
    tau = 20 * 25 * fscale[1][edgeR * dg_npf];
  }

  for(int m = 0; m < dg_np; m++) {
    for(int n = 0; n < dg_np; n++) {
      // op col-major
      int c_ind = m + n * dg_np;
      // op row-major
      // int c_ind = m * dg_np + n;
      op2L[c_ind] = 0.0;
      op2R[c_ind] = 0.0;
      for(int k = 0; k < dg_gf_np; k++) {
        // Dx' and Dy'
        int a_ind = m * dg_gf_np + k;
        // Dx and Dy
        int b_ind = n * dg_gf_np + k;

        int b_indP;
        if(*reverse) {
          b_indP = n * dg_gf_np + dg_gf_np - k - 1;
        } else {
          b_indP = n * dg_gf_np + k;
        }

        op1L[c_ind] += 0.5 * gaussW[k] * sJ[0][exIndL + k] * tau * gVML[a_ind] * gVML[b_ind];
        op1L[c_ind] += -0.5 * gaussW[k] * sJ[0][exIndL + k] * gVML[a_ind] * mDL[b_ind];
        op1L[c_ind] += -0.5 * gaussW[k] * sJ[0][exIndL + k] * mDL[a_ind] * gVML[b_ind];
        op1R[c_ind] += 0.5 * gaussW[k] * sJ[1][exIndR + k] * tau * gVMR[a_ind] * gVMR[b_ind];
        op1R[c_ind] += -0.5 * gaussW[k] * sJ[1][exIndR + k] * gVMR[a_ind] * mDR[b_ind];
        op1R[c_ind] += -0.5 * gaussW[k] * sJ[1][exIndR + k] * mDR[a_ind] * gVMR[b_ind];

        op2L[c_ind] += -0.5 * gaussW[k] * sJ[0][exIndL + k] * tau * gVML[a_ind] * gVMR[b_indP];
        op2L[c_ind] += -0.5 * gaussW[k] * sJ[0][exIndL + k] * gVML[a_ind] * pDL[b_ind];
        op2L[c_ind] += 0.5 * gaussW[k] * sJ[0][exIndL + k] * mDL[a_ind] * gVMR[b_indP];
        op2R[c_ind] += -0.5 * gaussW[k] * sJ[1][exIndR + k] * tau * gVMR[a_ind] * gVML[b_indP];
        op2R[c_ind] += -0.5 * gaussW[k] * sJ[1][exIndR + k] * gVMR[a_ind] * pDR[b_ind];
        op2R[c_ind] += 0.5 * gaussW[k] * sJ[1][exIndR + k] * mDR[a_ind] * gVML[b_indP];
      }
    }
  }
}
