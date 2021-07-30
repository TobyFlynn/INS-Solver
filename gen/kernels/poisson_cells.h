inline void poisson_cells(const double *u, const double *op, double *rhs) {
  for(int m = 0; m < 6; m++) {
    int ind = m * 6;
    rhs[m] = 0.0;
    for(int n = 0; n < 6; n++) {
      rhs[m] += op[ind + n] * u[n];
    }
  }
}
