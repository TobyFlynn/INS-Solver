inline void poisson_cells(const double *u, const double *op, double *rhs) {
  for(int m = 0; m < 10; m++) {
    int ind = m * 10;
    rhs[m] = 0.0;
    for(int n = 0; n < 10; n++) {
      rhs[m] += op[ind + n] * u[n];
    }
  }
}
