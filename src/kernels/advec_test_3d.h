inline void advec_test_3d(const double *x, const double *y, const double *z,
                          double *u, double *v, double *w, double *val) {
  for(int i = 0; i < DG_NP; i++) {
    u[i] = 1.0;
    v[i] = 1.0;
    w[i] = 1.0;
    val[i] = fmax(-(sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]) - 3.0), 0.0);
  }
}
