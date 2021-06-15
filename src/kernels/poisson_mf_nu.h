inline void poisson_mf_nu(const double *nu, double *dudx, double *dudy) {
  for(int i = 0; i < 15; i++) {
    dudx[i] = dudx[i] * nu[i];
    dudy[i] = dudy[i] * nu[i];
  }
}
