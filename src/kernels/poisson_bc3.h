inline void poisson_bc3(const double *bc, double *rhs) {
  for(int i = 0; i < 15; i++) {
    rhs[i] -= bc[i];
  }
}
