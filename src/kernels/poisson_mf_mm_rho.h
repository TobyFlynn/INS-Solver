inline void poisson_mf_mm_rho(const double *rho, const double *u, double *tmp) {
  for(int i = 0; i < 15; i++) {
    tmp[i] = u[i] * rho[i];
  }
}
