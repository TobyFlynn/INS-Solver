inline void poisson_mf2_bc(const double *tol, const int *bedge_type, const int *bedgeNum,
                           const int *d0, const int *d1, const int *d2,
                           const double *mD0, const double *mD1, const double *mD2,
                           const double *sJ, const double *tau, double *op) {
  const double *gVM, *mD;
  if(*bedgeNum == 0) {
    gVM = gFInterp0_g;
    mD  = mD0;
  } else if(*bedgeNum == 1) {
    gVM = gFInterp1_g;
    mD  = mD1;
  } else {
    gVM = gFInterp2_g;
    mD  = mD2;
  }

  if(*bedge_type == *d0 || *bedge_type == *d1 || *bedge_type == *d2) {
    for(int j = 0; j < DG_GF_NP * DG_NP; j++) {
      int indT = (j % DG_GF_NP) * DG_NP + (j / DG_GF_NP);
      int indT_col = j;
      int col  = j % DG_GF_NP;
      int row  = j / DG_GF_NP;
      double val = gaussW_g[j % DG_GF_NP] * sJ[*bedgeNum * DG_GF_NP + (j % DG_GF_NP)] * tau[*bedgeNum];
      val *= gVM[indT_col];
      val -= mD[indT] * gaussW_g[j % DG_GF_NP] * sJ[*bedgeNum * DG_GF_NP + (j % DG_GF_NP)];
      if(fabs(val) > *tol)
        op[row * DG_GF_NP + col] += val;
    }
  } else {
    for(int j = 0; j < DG_GF_NP * DG_NP; j++) {
      int indT = (j % DG_GF_NP) * DG_NP + (j / DG_GF_NP);
      int indT_col = j;
      int col  = j % DG_GF_NP;
      int row  = j / DG_GF_NP;
      double val = gaussW_g[j % DG_GF_NP] * sJ[*bedgeNum * DG_GF_NP + (j % DG_GF_NP)];
      val *= gVM[indT_col];
      if(fabs(val) > *tol)
        op[row * DG_GF_NP + col] += val;
    }
  }
}
