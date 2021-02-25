inline void viscosity_rhs(const double *factor, const double *J, double *vRHS0,
                          double *vRHS1) {

  for(int i = 0; i < 15; i++) {
    vRHS0[i] = *factor * J[i] * vRHS0[i];
    vRHS1[i] = *factor * J[i] * vRHS1[i];
  }
}
