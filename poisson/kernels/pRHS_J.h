inline void pRHS_J(const double *J, double *rhs) {
  for(int i = 0; i < 15; i++) {
    rhs[i] *= J[i];
  }
}
