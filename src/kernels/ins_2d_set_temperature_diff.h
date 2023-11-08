inline void ins_2d_set_temperature_diff(const DG_FP *rho, DG_FP *coeff) {
  for(int i = 0; i < DG_NP; i++) {
    coeff[i] = 1.0 / (peclet * rho[i]);
  }
}