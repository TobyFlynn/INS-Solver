inline void poisson_test_init(const double *x, const double *y, double *ex,
                              double *rhs, double *d, double *nx, double *ny) {
  for(int i = 0; i < 15; i++) {
    ex[i] = 0.0;
    double x1 = x[i];
    double y1 = y[i];
    // rhs[i] = -2.0 * (2.0 * (y[i] * y[i] * y[i]) - 3 * (y[i] * y[i]) + 1) + 6.0 * (1 - (x[i] * x[i])) * (2.0 * y[i] - 1.0);
    rhs[i] = 6 * x1 * y1 * (1.0 - y1) - 2.0 * x1 * x1 * x1;
    d[i] = 0.0;
    nx[i] = 0.0;
    ny[i] = 0.0;
  }
}
