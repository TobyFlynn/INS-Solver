inline void poisson_tmp(const double *op, double *tmp) {
  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind  = m * DG_NP + n;
      int indT = n * DG_NP + m;
      tmp[ind] = op[indT];
    }
  }
}
