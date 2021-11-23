inline void poisson_cells(const double *u, const double *op, double *rhs) {
  for(int m = 0; m < DG_NP; m++) {
    rhs[m] = 0.0;
    for(int n = 0; n < DG_NP; n++) {
      int ind = m + n * DG_NP;
      rhs[m] += op[ind] * u[n];
    }
  }
}
