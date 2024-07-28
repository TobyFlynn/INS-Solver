inline void measure_mass_of_phases(const DG_FP *alpha, const DG_FP *s, const DG_FP *rho,
                                   DG_FP *vol_0, DG_FP *vol_1, DG_FP *mass_0, 
                                   DG_FP *mass_1, DG_FP *mass_total) {
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    DG_FP step = tanh(PI * s[i] / *alpha);
    vol_0[i] = fmax(step, 0.0);
    vol_1[i] = fmin(step, 0.0);
    mass_0[i] = step >= 0.0 ? rho[i] : 0.0;
    mass_1[i] = step < 0.0 ? rho[i] : 0.0;
    mass_total[i] = rho[i];
  }
}