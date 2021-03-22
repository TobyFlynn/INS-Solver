inline void viscosity_rhs(const double *factor, const double *J, double *vRHS0,
                          double *vRHS1, double *bcx, double *bcy) {

  for(int i = 0; i < 15; i++) {
    // vRHS0[i] = *factor * J[i] * vRHS0[i];
    // vRHS1[i] = *factor * J[i] * vRHS1[i];
    vRHS0[i] = (*factor) * vRHS0[i];
    vRHS1[i] = (*factor) * vRHS1[i];
  }

  // for(int i = 0; i < 21; i++) {
  //   bcx[i] *= -1.0;
  //   bcy[i] *= -1.0;
  // }
}
