inline void poisson_mf2_opf(const double *tol, const int *edgeNum, const double *gop0L,
                            const double *gop1L, const double *gop2L, const double *gopf0L,
                            const double *gopf1L, const double *gopf2L, double *op20L,
                            double *op21L, double *op22L, double *op1L,
                            const double *gop0R, const double *gop1R, const double *gop2R,
                            const double *gopf0R, const double *gopf1R, const double *gopf2R,
                            double *op20R, double *op21R, double *op22R, double *op1R) {
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  if(edgeL == 0) {
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gop0L[colInd];
        // if(fabs(val) > *tol) {
          op1L[ind] += val;
        // }
        val = -0.5 * gopf0L[colInd];
        // if(fabs(val) > *tol) {
          op20L[ind] += val;
        // }
      }
    }
  } else if(edgeL == 1) {
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gop1L[colInd];
        // if(fabs(val) > *tol) {
          op1L[ind] += val;
        // }
        val = -0.5 * gopf1L[colInd];
        // if(fabs(val) > *tol) {
          op21L[ind] += val;
        // }
      }
    }
  } else {
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gop2L[colInd];
        // if(fabs(val) > *tol) {
          op1L[ind] += val;
        // }
        val = -0.5 * gopf2L[colInd];
        // if(fabs(val) > *tol) {
          op22L[ind] += val;
        // }
      }
    }
  }

  if(edgeR == 0) {
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gop0R[colInd];
        // if(fabs(val) > *tol) {
          op1R[ind] += val;
        // }
        val = -0.5 * gopf0R[colInd];
        // if(fabs(val) > *tol) {
          op20R[ind] += val;
        // }
      }
    }
  } else if(edgeR == 1) {
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gop1R[colInd];
        // if(fabs(val) > *tol) {
          op1R[ind] += val;
        // }
        val = -0.5 * gopf1R[colInd];
        // if(fabs(val) > *tol) {
          op21R[ind] += val;
        // }
      }
    }
  } else {
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gop2R[colInd];
        // if(fabs(val) > *tol) {
          op1R[ind] += val;
        // }
        val = -0.5 * gopf2R[colInd];
        // if(fabs(val) > *tol) {
          op22R[ind] += val;
        // }
      }
    }
  }
}
