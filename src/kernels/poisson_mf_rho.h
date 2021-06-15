inline void poisson_mf_rho(const double *rho, double *dudx, double *dudy) {
  for(int i = 0; i < 15; i++) {
    dudx[i] = dudx[i] / rho[i];
    dudy[i] = dudy[i] / rho[i];
  }
}
