inline void init_nu_rho(const int *p, double *nu, double *rho) {
  // Get constants for this element's order
  const int dg_np  = DG_CONSTANTS[(*p - 1) * 5];
  
  for(int i = 0; i < dg_np; i++) {
    nu[i] = nu0;
    rho[i] = rh0;
  }
}
