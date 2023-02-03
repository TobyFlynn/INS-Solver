inline void advec_test_ic(const double *x, const double *y, double *val,
                          double *u, double *v) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    val[i] = sqrt(x[i] * x[i] + y[i] * y[i]) - 2.0;
    val[i] = fmax(-1.0 * val[i], 0.0);
    u[i]   = 1.0;
    v[i]   = 1.0;
  }
}
