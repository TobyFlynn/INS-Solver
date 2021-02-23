inline void poisson_test_set_rhs(const double *J, const double *ex, double *rhs) {
/*
  for(int i = 0; i < 15; i++) {
    if(ex[i] < 0.0)
      rhs[FMASK[i]] = 0.0;
  }
*/
  for(int i = 0; i < 15; i++) {
    rhs[i] *= J[i];
  }
}
