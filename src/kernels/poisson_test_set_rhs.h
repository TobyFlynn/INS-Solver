inline void poisson_test_set_rhs(const double *ex, double *rhs) {
  for(int i = 0; i < 15; i++) {
    if(ex[i] < 0.0)
      rhs[FMASK[i]] = 0.0;
  }
}
