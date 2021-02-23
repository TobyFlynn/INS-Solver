inline void viscosity_rhs(const double *factor, const double *J, const double *qtt0,
                          const double *qtt1, double *vRHS0, double *vRHS1) {

  for(int i = 0; i < 15; i++) {
    vRHS0[i] = *factor * J[i] * qtt0[i];
    vRHS1[i] = *factor * J[i] * qtt1[i];
    // vRHS0[i] = *factor * qtt0[i];
    // vRHS1[i] = *factor * qtt1[i];
  }
}
