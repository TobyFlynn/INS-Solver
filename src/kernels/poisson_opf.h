inline void poisson_opf(const double *tol, const int *edgeNum, const double *gop0L,
                            const double *gop1L, const double *gop2L, const double *gopf0L,
                            const double *gopf1L, const double *gopf2L, double *op2L,
                            double *op1L,
                            const double *gop0R, const double *gop1R, const double *gop2R,
                            const double *gopf0R, const double *gopf1R, const double *gopf2R,
                            double *op2R, double *op1R) {
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  const double *gopL, *gopR, *gopFL, *gopFR;
  if(edgeL == 0) {
    gopL  = gop0L;
    gopFL = gopf0L;
  } else if(edgeL == 1) {
    gopL  = gop1L;
    gopFL = gopf1L;
  } else {
    gopL  = gop2L;
    gopFL = gopf2L;
  }

  if(edgeR == 0) {
    gopR  = gop0R;
    gopFR = gopf0R;
  } else if(edgeR == 1) {
    gopR  = gop1R;
    gopFR = gopf1R;
  } else {
    gopR  = gop2R;
    gopFR = gopf2R;
  }

  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind = m * DG_NP + n;
      int colInd = n * DG_NP + m;

      double val = 0.5 * gopL[colInd];
      if(fabs(val) > *tol)
        op1L[ind] += val;
      val = -0.5 * gopFL[colInd];
      op2L[ind] = val;

      val = 0.5 * gopR[colInd];
      if(fabs(val) > *tol)
        op1R[ind] += val;
      val = -0.5 * gopFR[colInd];
      op2R[ind] = val;
    }
  }
}
