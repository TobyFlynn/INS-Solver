inline void mp_ins_pressure_fact_2d(const DG_FP *rho, DG_FP *fact) {
  for(int i = 0; i < DG_NP; i++) {
    fact[i] = 1.0 / rho[i];
  }
}
