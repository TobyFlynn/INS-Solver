inline void init_nu_rho(double *nu, double *rho) {
  for(int i = 0; i < DG_NP; i++) {
    nu[i] = nu0;
    rho[i] = 1.0;
  }
}
