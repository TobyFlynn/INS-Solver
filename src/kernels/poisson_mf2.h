inline void poisson_mf2(const double *u, const double *op, double *rhs) {
  for(int m = 0; m < 15; m++) {
    int ind = m * 15;
    double val = 0.0;
    for(int n = 0; n < 15; n++) {
      val += op[ind + n] * u[n];
    }
    rhs[m] = val;
  }
}
