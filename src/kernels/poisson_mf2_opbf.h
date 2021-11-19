inline void poisson_mf2_opbf(const double *tol, const int *btype, const int *edgeNum,
                             const int *d0, const int *d1, const int *d2, const double *gop0,
                             const double *gop1, const double *gop2, double *op1) {
  if(*btype == *d0 || *btype == *d1 || *btype == *d2) {
    const double *gop;
    if(*edgeNum == 0) {
      gop = gop0;
    } else if(*edgeNum == 1) {
      gop = gop1;
    } else {
      gop = gop2;
    }

    for(int m = 0; m < DG_NP; m++) {
      for(int n = 0; n < DG_NP; n++) {
        int ind = m * DG_NP + n;
        int colInd = n * DG_NP + m;
        double val = gop[colInd];
        if(fabs(val) > *tol) {
          op1[ind] += val;
        }
      }
    }
  }
}
