inline void poisson_cells(const double *u, const double *op, double *rhs) {
  for(int m = 0; m < 15; m++) {
    int ind = m * 15;
    rhs[m] = 0.0;
    for(int n = 0; n < 15; n++) {
      rhs[m] += op[ind + n] * u[n];
    }
  }
}
