inline void pRHS_du(const double *U, const double *exU, double *du) {
  for(int i = 0; i < 15; i++) {
    du[i] = U[i] - exU[i];
  }
}
