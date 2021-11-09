inline void ls_fluid_prop(const double *step, double *nu, double *rho) {
  for(int i = 0; i < DG_NP; i++) {
    nu[i] = 0.5 * nu0 * (1.0 + step[i]) + 0.5 * nu1 * (1.0 - step[i]);
    rho[i] = 0.5 * rho0 * (1.0 + step[i]) + 0.5 * rho1 * (1.0 - step[i]);
  }
}
