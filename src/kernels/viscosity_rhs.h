inline void viscosity_rhs(const double *factor, const double *rho,
                          const double *qtt0, const double *qtt1,
                          double *vRHS0, double *vRHS1) {
  for(int i = 0; i < DG_NP; i++) {
    vRHS0[i] = (*factor) * rho[i] * qtt0[i];
    vRHS1[i] = (*factor) * rho[i] * qtt1[i];
  }
}
