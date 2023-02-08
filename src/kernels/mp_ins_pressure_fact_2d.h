inline void mp_ins_pressure_fact_2d(const double *rho, double *fact) {
  for(int i = 0; i < DG_NP; i++) {
    fact[i] = 1.0 / rho[i];
  }
}
