inline void viscosity_rhs(const double *factor, const double *nu, double *vRHS0,
                          double *vRHS1, double *vFactor) {
  for(int i = 0; i < 15; i++) {
    vFactor[i] = (*factor) / nu[i];
    vRHS0[i] = vFactor[i] * vRHS0[i];
    vRHS1[i] = vFactor[i] * vRHS1[i];
  }
}
