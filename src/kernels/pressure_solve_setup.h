inline void pressure_solve_setup(const double *rho, double *factor) {
  for(int i = 0; i < 15; i++) {
    factor[i] = 1.0 / rho[i];
  }
}
