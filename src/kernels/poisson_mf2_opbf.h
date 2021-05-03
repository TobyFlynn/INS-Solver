inline void poisson_mf2_opbf(const double *tol, const int *btype, const int *edgeNum,
                             const int *d0, const int *d1, const int *d2, const double *gop0,
                             const double *gop1, const double *gop2, double *op1) {
  if(*btype == *d0 || *btype == *d1 || *btype == *d2) {
    if(*edgeNum == 0) {
      for(int m = 0; m < 15; m++) {
        for(int n = 0; n < 15; n++) {
          int ind = m * 15 + n;
          int colInd = n * 15 + m;
          double val = gop0[colInd];
          if(fabs(val) > *tol) {
            op1[ind] += val;
          }
        }
      }
    } else if(*edgeNum == 1) {
      for(int m = 0; m < 15; m++) {
        for(int n = 0; n < 15; n++) {
          int ind = m * 15 + n;
          int colInd = n * 15 + m;
          double val = gop1[colInd];
          if(fabs(val) > *tol) {
            op1[ind] += val;
          }
        }
      }
    } else {
      for(int m = 0; m < 15; m++) {
        for(int n = 0; n < 15; n++) {
          int ind = m * 15 + n;
          int colInd = n * 15 + m;
          double val = gop2[colInd];
          if(fabs(val) > *tol) {
            op1[ind] += val;
          }
        }
      }
    }
  }
}
