inline void poisson_gauss_grad_b(const int *btype, const int *edgeNum,
                                 const int *d0, const int *d1, const int *d2,
                                 const double *x, const double *y,
                                 const double *sJ, const double *nx,
                                 const double *ny, const double *fscale,
                                 double *op1, double *op_bc) {
  const double *gVM;
  if(*edgeNum == 0) {
    gVM = gFInterp0_g;
  } else if(*edgeNum == 1) {
    gVM = gFInterp1_g;
  } else {
    gVM = gFInterp2_g;
  }

  if(*d0 == *btype || *d1 == *btype || *d2 == *btype) {
    const double *gDr, *gDs, *gVM;
    if(*edgeNum == 0) {
      gDr = gF0Dr_g;
      gDs = gF0Ds_g;
      gVM = gFInterp0_g;
    } else if(*edgeNum == 1) {
      gDr = gF1Dr_g;
      gDs = gF1Ds_g;
      gVM = gFInterp1_g;
    } else {
      gDr = gF2Dr_g;
      gDs = gF2Ds_g;
      gVM = gFInterp2_g;
    }

    double rx[DG_GF_NP], sx[DG_GF_NP], ry[DG_GF_NP], sy[DG_GF_NP];

    for(int m = 0; m < DG_GF_NP; m++) {
      rx[m] = 0.0;
      sx[m] = 0.0;
      ry[m] = 0.0;
      sy[m] = 0.0;
      for(int n = 0; n < DG_NP; n++) {
        int ind = m + n * DG_GF_NP;
        rx[m] += gDr[ind] * x[n];
        sx[m] += gDs[ind] * x[n];
        ry[m] += gDr[ind] * y[n];
        sy[m] += gDs[ind] * y[n];
      }
      double J = -sx[m] * ry[m] + rx[m] * sy[m];
      double rx_n = sy[m] / J;
      double sx_n = -ry[m] / J;
      double ry_n = -sx[m] / J;
      double sy_n = rx[m] / J;
      rx[m] = rx_n;
      sx[m] = sx_n;
      ry[m] = ry_n;
      sy[m] = sy_n;
    }

    const int exInd = *edgeNum * DG_GF_NP;
    double mD[DG_GF_NP * DG_NP];
    for(int m = 0; m < DG_GF_NP; m++) {
      for(int n = 0; n < DG_NP; n++) {
        int ind = m + n * DG_GF_NP;

        double Dx = rx[m] * gDr[ind] + sx[m] * gDs[ind];
        double Dy = ry[m] * gDr[ind] + sy[m] * gDs[ind];
        mD[ind]   = nx[exInd + m] * Dx + ny[exInd + m] * Dy;
      }
    }

    double tau = 20 * 25 * fscale[*edgeNum * DG_NPF];
    // Main matrix
    for(int m = 0; m < DG_NP; m++) {
      for(int n = 0; n < DG_NP; n++) {
        // op col-major
        // int c_ind = m + n * DG_NP;
        // op row-major
        int c_ind = m * DG_NP + n;
        for(int k = 0; k < DG_GF_NP; k++) {
          // Dx' and Dy'
          int a_ind = m * DG_GF_NP + k;
          // Dx and Dy
          int b_ind = n * DG_GF_NP + k;
          op1[c_ind] += gaussW_g[k] * sJ[exInd + k] * tau * gVM[a_ind] * gVM[b_ind];
          op1[c_ind] += -gaussW_g[k] * sJ[exInd + k] * gVM[a_ind] * mD[b_ind];
          op1[c_ind] += -gaussW_g[k] * sJ[exInd + k] * mD[a_ind] * gVM[b_ind];
        }
      }
    }

    // Apply BC matrix
    for(int j = 0; j < DG_GF_NP * DG_NP; j++) {
      int indT_col = j;
      int col  = j % DG_GF_NP;
      int row  = j / DG_GF_NP;
      double val = gaussW_g[j % DG_GF_NP] * sJ[*edgeNum * DG_GF_NP + (j % DG_GF_NP)] * tau;
      val *= gVM[indT_col];
      val -= mD[indT_col] * gaussW_g[j % DG_GF_NP] * sJ[*edgeNum * DG_GF_NP + (j % DG_GF_NP)];
      op_bc[row * DG_GF_NP + col] = val;
    }
  } else {
    // Nothing for main matrix
    // Apply BC matrix
    for(int j = 0; j < DG_GF_NP * DG_NP; j++) {
      int indT_col = j;
      int col  = j % DG_GF_NP;
      int row  = j / DG_GF_NP;
      double val = gaussW_g[j % DG_GF_NP] * sJ[*edgeNum * DG_GF_NP + (j % DG_GF_NP)];
      val *= gVM[indT_col];
      op_bc[row * DG_GF_NP + col] = val;
    }
  }
}
