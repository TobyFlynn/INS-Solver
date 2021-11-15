inline void advection_flux(const int *, const double *u, const double *v,
                           double *f0, double *f1, double *f2, double *f3) {
  // Get constants for this element's order
  const int dg_np  = DG_CONSTANTS[(*p - 1) * 5];

  for(int i = 0; i < dg_np; i++) {
    f0[i] = u[i] * u[i];
    f1[i] = u[i] * v[i];
    f2[i] = u[i] * v[i];
    f3[i] = v[i] * v[i];
  }
}
