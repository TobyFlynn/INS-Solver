inline void poisson_cells(const double *u, const double *op, double *rhs) {
  for(int m = 0; m < 3; m++) {
    int ind = m * 3;
    rhs[m] = 0.0;
    for(int n = 0; n < 3; n++) {
      rhs[m] += op[ind + n] * u[n];
    }
  }
}
