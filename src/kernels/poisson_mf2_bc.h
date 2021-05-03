inline void poisson_mf2_bc(const double *tol, const int *bedge_type, const int *bedgeNum,
                           const int *d0, const int *d1, const int *d2,
                           const double *mD0, const double *mD1, const double *mD2,
                           const double *sJ, const double *tau, double *op) {
  if(*bedge_type == *d0 || *bedge_type == *d1 || *bedge_type == *d2) {
    for(int j = 0; j < 7 * 15; j++) {
      int indT = (j % 7) * 15 + (j / 7);
      int col  = j % 7;
      int row  = j / 7;
      double val = gaussW_g[j % 7] * sJ[*bedgeNum * 7 + (j % 7)] * tau[*bedgeNum];
      double mD;
      if(*bedgeNum == 0) {
        val *= gFInterp0_g[indT];
        mD = mD0[indT];
      } else if(*bedgeNum == 1) {
        val *= gFInterp1_g[indT];
        mD = mD1[indT];
      } else {
        val *= gFInterp2_g[indT];
        mD = mD2[indT];
      }
      val -= mD * gaussW_g[j % 7] * sJ[*bedgeNum * 7 + (j % 7)];
      if(fabs(val) > *tol)
        op[row * 7 + col] += val;
    }
  } else {
    for(int j = 0; j < 7 * 15; j++) {
      int indT = (j % 7) * 15 + (j / 7);
      int col  = j % 7;
      int row  = j / 7;
      double val = gaussW_g[j % 7] * sJ[*bedgeNum * 7 + (j % 7)];
      if(*bedgeNum == 0) {
        val *= gFInterp0_g[indT];
      } else if(*bedgeNum == 1) {
        val *= gFInterp1_g[indT];
      } else {
        val *= gFInterp2_g[indT];
      }
      if(fabs(val) > *tol)
        op[row * 7 + col] += val;
    }
  }
}
