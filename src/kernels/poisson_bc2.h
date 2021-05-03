inline void poisson_bc2(const double *sJ, const double *tau, const double *bc,
                        double *gtau) {
  for(int i = 0; i < 21; i++) {
    gtau[i] = tau[i / 7] * gaussW_g[i % 7] * sJ[i] * bc[i];
  }
}
