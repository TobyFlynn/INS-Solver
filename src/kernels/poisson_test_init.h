inline void poisson_test_init(double *ex, double *rhs) {
  for(int i = 0; i < 15; i++) {
    ex[i] = 0.0;
    rhs[i] = 1.0;
  }
}
