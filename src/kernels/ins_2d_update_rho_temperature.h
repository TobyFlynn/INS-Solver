inline void ins_2d_update_rho_temperature(const DG_FP *temperature, DG_FP *rho) {
  for(int i = 0; i < DG_NP; i++) {
    // rho[i] = rho_0 - coeff_thermal_expan * rho_0 * (temperature[i] - temperature_0);
    // Non dimensionalise to the following
    rho[i] = 1.0 - coeff_thermal_expan * temperature[i];
  }
}