inline void poisson_test_tmp(double *rhs) {
  for(int i = 0; i < 15; i++) {
    rhs[i] *= -1.0;
  }
}
