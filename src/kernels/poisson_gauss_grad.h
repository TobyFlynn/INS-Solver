inline void poisson_gauss_grad(const int *edgeNum, const bool *reverse,
                               const double **x, const double **y,
                               const double **sJ, const double **nx,
                               const double **ny, const double **fscale,
                               double *op1L, double *op1R, double *op2L,
                               double *op2R) {
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  const double *gDrL, *gDsL, *gDrR, *gDsR, *gVML, *gVMR;

  if(edgeL == 0) {
    gDrL = gF0Dr_g;
    gDsL = gF0Ds_g;
    gVML = gFInterp0_g;
  } else if(edgeL == 1) {
    gDrL = gF1Dr_g;
    gDsL = gF1Ds_g;
    gVML = gFInterp1_g;
  } else {
    gDrL = gF2Dr_g;
    gDsL = gF2Ds_g;
    gVML = gFInterp2_g;
  }

  if(edgeR == 0) {
    gDrR = gF0Dr_g;
    gDsR = gF0Ds_g;
    gVMR = gFInterp0_g;
  } else if(edgeR == 1) {
    gDrR = gF1Dr_g;
    gDsR = gF1Ds_g;
    gVMR = gFInterp1_g;
  } else {
    gDrR = gF2Dr_g;
    gDsR = gF2Ds_g;
    gVMR = gFInterp2_g;
  }

  double rxL[DG_GF_NP], sxL[DG_GF_NP], ryL[DG_GF_NP], syL[DG_GF_NP];
  double rxR[DG_GF_NP], sxR[DG_GF_NP], ryR[DG_GF_NP], syR[DG_GF_NP];

  for(int m = 0; m < DG_GF_NP; m++) {
    rxL[m] = 0.0;
    sxL[m] = 0.0;
    ryL[m] = 0.0;
    syL[m] = 0.0;
    rxR[m] = 0.0;
    sxR[m] = 0.0;
    ryR[m] = 0.0;
    syR[m] = 0.0;
    for(int n = 0; n < DG_NP; n++) {
      int ind = m + n * DG_GF_NP;
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

  const int exIndL = edgeL * DG_GF_NP;
  const int exIndR = edgeR * DG_GF_NP;
  double mDL[DG_GF_NP * DG_NP], mDR[DG_GF_NP * DG_NP], pDL[DG_GF_NP * DG_NP], pDR[DG_GF_NP * DG_NP];
  for(int m = 0; m < DG_GF_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind = m + n * DG_GF_NP;
      int p_ind, p_norm_indL, p_norm_indR;
      if(*reverse) {
        p_ind = DG_GF_NP - 1 - m + n * DG_GF_NP;
        p_norm_indL = exIndL + DG_GF_NP - 1 - m;
        p_norm_indR = exIndR + DG_GF_NP - 1 - m;
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
  if(fscale[0][edgeL * DG_NPF] > fscale[1][edgeR * DG_NPF]) {
    tau = 20 * 25 * fscale[0][edgeL * DG_NPF];
  } else {
    tau = 20 * 25 * fscale[1][edgeR * DG_NPF];
  }

  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      // op col-major
      int c_ind = m + n * DG_NP;
      // op row-major
      // int c_ind = m * DG_NP + n;
      op2L[c_ind] = 0.0;
      op2R[c_ind] = 0.0;
      for(int k = 0; k < DG_GF_NP; k++) {
        // Dx' and Dy'
        int a_ind = m * DG_GF_NP + k;
        // Dx and Dy
        int b_ind = n * DG_GF_NP + k;

        int b_indP;
        if(*reverse) {
          b_indP = n * DG_GF_NP + DG_GF_NP - k - 1;
        } else {
          b_indP = n * DG_GF_NP + k;
        }

        op1L[c_ind] += 0.5 * gaussW_g[k] * sJ[0][exIndL + k] * tau * gVML[a_ind] * gVML[b_ind];
        op1L[c_ind] += -0.5 * gaussW_g[k] * sJ[0][exIndL + k] * gVML[a_ind] * mDL[b_ind];
        op1L[c_ind] += -0.5 * gaussW_g[k] * sJ[0][exIndL + k] * mDL[a_ind] * gVML[b_ind];
        op1R[c_ind] += 0.5 * gaussW_g[k] * sJ[1][exIndR + k] * tau * gVMR[a_ind] * gVMR[b_ind];
        op1R[c_ind] += -0.5 * gaussW_g[k] * sJ[1][exIndR + k] * gVMR[a_ind] * mDR[b_ind];
        op1R[c_ind] += -0.5 * gaussW_g[k] * sJ[1][exIndR + k] * mDR[a_ind] * gVMR[b_ind];

        op2L[c_ind] += -0.5 * gaussW_g[k] * sJ[0][exIndL + k] * tau * gVML[a_ind] * gVMR[b_indP];
        op2L[c_ind] += -0.5 * gaussW_g[k] * sJ[0][exIndL + k] * gVML[a_ind] * pDL[b_ind];
        op2L[c_ind] += 0.5 * gaussW_g[k] * sJ[0][exIndL + k] * mDL[a_ind] * gVMR[b_indP];
        op2R[c_ind] += -0.5 * gaussW_g[k] * sJ[1][exIndR + k] * tau * gVMR[a_ind] * gVML[b_indP];
        op2R[c_ind] += -0.5 * gaussW_g[k] * sJ[1][exIndR + k] * gVMR[a_ind] * pDR[b_ind];
        op2R[c_ind] += 0.5 * gaussW_g[k] * sJ[1][exIndR + k] * mDR[a_ind] * gVML[b_indP];
      }
    }
  }
}
