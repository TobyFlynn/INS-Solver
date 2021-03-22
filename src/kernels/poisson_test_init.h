inline void poisson_test_init(const double *x, const double *y, double *ex,
                              double *rhs) {
  for(int i = 0; i < 15; i++) {
    ex[i] = 0.0;
    double x1 = x[i];
    double y1 = y[i];
    rhs[i] = -2.0 * (2.0 * (y1 * y1 * y1) - 3 * (y1 * y1) + 1) + 6.0 * (1 - (x1 * x1)) * (2.0 * y1 - 1.0);
    // rhs[i] = 6 * x1 * y1 * (1.0 - y1) - 2.0 * x1 * x1 * x1;
    // rhs[i] = 1.0;
  }
}
