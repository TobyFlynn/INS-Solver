inline void ls_test_ic(const double *x, const double *y,
                          double *u, double *v) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    u[i]   = 1.0;
    v[i]   = 1.0;
  }
}
