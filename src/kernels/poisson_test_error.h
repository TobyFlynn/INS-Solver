inline void poisson_test_error(const double *x, const double *y,
                               const double *sol, double *err) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < 15; i++) {
    double exact = - (8.0/(PI*PI)) * (sin(x[i])*sin(y[i]) + (1.0/15.0)*sin(x[i])*sin(3.0*y[i]) + (1.0/15.0)*sin(3.0*x[i])*sin(y[i]));
    err[i] = fabs(sol[i] - exact);
  }
}
