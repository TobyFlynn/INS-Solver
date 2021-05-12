inline void poisson_test_error(const double *x, const double *y,
                               const double *sol, double *err) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < 15; i++) {
    double x1 = x[i];
    double y1 = y[i];
    // double exact = - (8.0/(PI*PI)) * (sin(x[i])*sin(y[i]) + (1.0/15.0)*sin(x[i])*sin(3.0*y[i]) + (1.0/15.0)*sin(3.0*x[i])*sin(y[i]) + (1.0/81.0)*sin(3.0*x[i])*sin(3.0*y[i]));
    // exact += - (8.0/(PI*PI)) * ((1.0/65.0)*sin(5.0*x[i])*sin(y[i]) + (1.0/65.0)*sin(x[i])*sin(5.0*y[i]) + (1.0/255.0)*sin(5.0*x[i])*sin(3.0*y[i]) + (1.0/255.0)*sin(3.0*x[i])*sin(5.0*y[i]) + (1.0/625.0)*sin(5.0*x[i])*sin(5.0*y[i]));
    // double exact = (1.0 - (x[i] * x[i])) * (2.0 * (y[i] * y[i] * y[i]) - 3.0 * (y[i] * y[i]) + 1.0);
    double exact = y1 * (1.0 - y1) * x1 * x1 * x1;
    err[i] = fabs(sol[i] - exact);
    // err[i] = fabs(sol[i] - sol2[i]);
  }
}
