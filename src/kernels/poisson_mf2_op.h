inline void poisson_mf2_op(const double *cub_op, const double *tol, double *op1) {
  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind = m * DG_NP + n;
      int colInd = n * DG_NP + m;
      if(fabs(cub_op[colInd]) > *tol) {
        op1[ind] = cub_op[colInd];
      }
    }
  }
}
