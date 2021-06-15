inline void poisson_mf_bc(const int *bedgeType, const int *bedgeNum,
                          const int *d0, const int *d1, const int *d2,
                          const double *mD0, const double *mD1,
                          const double *mD2, const double *sJ,
                          const double *tau, const double *bc, double *rhs) {
  int exInd = 0;
  if(*bedgeNum == 1) exInd = 7;
  else if(*bedgeNum == 2) exInd = 2 * 7;

  if(*bedgeType == *d0 || *bedgeType == *d1 || *bedgeType == *d2) {
    // Temp rough work just to see if it works
    double op[15 * 7];
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
      // if(fabs(val) > 1e-15)
        op[row * 7 + col] += val;
    }

    for(int m = 0; m < 15; m++) {
      int ind = m * 7;
      double val = 0.0;
      for(int n = 0; n < 7; n++) {
        val += op[ind + n] * bc[exInd + n];
      }
      rhs[m] += val;
    }
  } else {
    double op[15 * 7];
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
      // if(fabs(val) > 1e-15)
        op[row * 7 + col] += val;
    }

    for(int m = 0; m < 15; m++) {
      int ind = m * 7;
      double val = 0.0;
      for(int n = 0; n < 7; n++) {
        val += op[ind + n] * bc[exInd + n];
      }
      rhs[m] += val;
    }
  }
}
