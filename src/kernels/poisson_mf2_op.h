inline void poisson_mf2_op(const double *cub_op, const double *tol, double *op1) {
  for(int m = 0; m < 15; m++) {
    for(int n = 0; n < 15; n++) {
      int ind = m * 15 + n;
      int colInd = n * 15 + m;
      if(fabs(cub_op[colInd]) > *tol) {
        op1[ind] = cub_op[colInd];
      }
    }
  }
}
