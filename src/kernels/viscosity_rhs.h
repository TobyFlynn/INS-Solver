inline void viscosity_rhs(const double *factor, const double *rho,
                          double *vRHS0, double *vRHS1) {
  for(int i = 0; i < DG_NP; i++) {
    vRHS0[i] = (*factor) * rho[i] * vRHS0[i];
    vRHS1[i] = (*factor) * rho[i] * vRHS1[i];
  }
}
