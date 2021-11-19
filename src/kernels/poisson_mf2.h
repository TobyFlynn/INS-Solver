inline void poisson_mf2(const double *u, const double *op, double *rhs) {
  for(int m = 0; m < DG_NP; m++) {
    int ind = m * DG_NP;
    double val = 0.0;
    for(int n = 0; n < DG_NP; n++) {
      val += op[ind + n] * u[n];
    }
    rhs[m] = val;
  }
}
