inline void poisson_mf_bc1(const double *tmp, double *b) {
  for(int i = 0; i < 15; i++) {
    b[i] += tmp[i];
  }
}
