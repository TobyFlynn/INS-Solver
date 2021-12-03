inline void poisson_gauss_grad_b(const int *p, const double *gF0Dr,
                                 const double *gF0Ds, const double *gF1Dr,
                                 const double *gF1Ds, const double *gF2Dr,
                                 const double *gF2Ds, const double *gFInterp0,
                                 const double *gFInterp1, const double *gFInterp2,
                                 const int *btype, const int *edgeNum,
                                 const int *d0, const int *d1, const int *d2,
                                 const double *x, const double *y,
                                 const double *sJ, const double *nx,
                                 const double *ny, const double *h,
                                 const double *factor, double *op1,
                                 double *op_bc) {
  const double *gVM;
  if(*edgeNum == 0) {
    gVM = &gFInterp0[(*p - 1) * DG_GF_NP * DG_NP];
  } else if(*edgeNum == 1) {
    gVM = &gFInterp1[(*p - 1) * DG_GF_NP * DG_NP];
  } else {
    gVM = &gFInterp2[(*p - 1) * DG_GF_NP * DG_NP];
  }

  // Get constants
  const int dg_np      = DG_CONSTANTS[(*p - 1) * 5];
  const int dg_npf     = DG_CONSTANTS[(*p - 1) * 5 + 1];
  const int dg_gf_np   = DG_CONSTANTS[(*p - 1) * 5 + 4];
  const double *gaussW = &gaussW_g[(*p - 1) * DG_GF_NP];

  if(*d0 == *btype || *d1 == *btype || *d2 == *btype) {
    const double *gDr, *gDs;
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

    double rx[DG_GF_NP], sx[DG_GF_NP], ry[DG_GF_NP], sy[DG_GF_NP];

    for(int m = 0; m < dg_gf_np; m++) {
      rx[m] = 0.0;
      sx[m] = 0.0;
      ry[m] = 0.0;
      sy[m] = 0.0;
      for(int n = 0; n < dg_np; n++) {
        int ind = m + n * dg_gf_np;
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

    const int exInd = *edgeNum * dg_gf_np;
    double mD[DG_GF_NP * DG_NP];
    for(int m = 0; m < dg_gf_np; m++) {
      for(int n = 0; n < dg_np; n++) {
        int ind = m + n * dg_gf_np;

        double Dx = rx[m] * gDr[ind] + sx[m] * gDs[ind];
        double Dy = ry[m] * gDr[ind] + sy[m] * gDs[ind];
        mD[ind]   = factor[exInd + m] * (nx[exInd + m] * Dx + ny[exInd + m] * Dy);
      }
    }

    double tau[DG_GF_NP];
    for(int i = 0; i < DG_GF_NP; i++) {
      int ind = *edgeNum * DG_GF_NP + i;
      tau[i] = (DG_ORDER + 1) * (DG_ORDER + 2) * (*h) * factor[ind];
    }

    // Main matrix
    for(int m = 0; m < dg_np; m++) {
      for(int n = 0; n < dg_np; n++) {
        // op col-major
        int c_ind = m + n * dg_np;
        // op row-major
        // int c_ind = m * dg_np + n;
        for(int k = 0; k < dg_gf_np; k++) {
          // Dx' and Dy'
          int a_ind = m * dg_gf_np + k;
          // Dx and Dy
          int b_ind = n * dg_gf_np + k;
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
      double val = gaussW[j % dg_gf_np] * sJ[*edgeNum * dg_gf_np + (j % dg_gf_np)] * tau[j % dg_gf_np];
      val *= gVM[indT_col];
      val -= mD[indT_col] * gaussW[j % dg_gf_np] * sJ[*edgeNum * dg_gf_np + (j % dg_gf_np)];
      op_bc[row + col * dg_np] = val;
    }
  } else {
    // Nothing for main matrix
    // Apply BC matrix
    for(int j = 0; j < dg_gf_np * dg_np; j++) {
      int indT_col = j;
      int col  = j % dg_gf_np;
      int row  = j / dg_gf_np;
      double val = gaussW[j % dg_gf_np] * sJ[*edgeNum * dg_gf_np + (j % dg_gf_np)];
      val *= gVM[indT_col];
      op_bc[row + col * dg_np] = val;
    }
  }
}
