inline void poisson_bc_J(const double *J, double *gx, double *gy) {
  for(int i = 0; i < 15; i++) {
    gx[i] = gx[i] / J[i];
    gy[i] = gy[i] / J[i];
  }
}
