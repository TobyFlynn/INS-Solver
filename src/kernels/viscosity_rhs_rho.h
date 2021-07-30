inline void viscosity_rhs_rho(const double *rho, double *vRHS0, double *vRHS1) {
  for(int i = 0; i < DG_NP; i++) {
    vRHS0[i] = rho[i] * vRHS0[i];
    vRHS1[i] = rho[i] * vRHS1[i];
  }
}
