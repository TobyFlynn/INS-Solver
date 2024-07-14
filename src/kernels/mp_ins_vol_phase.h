inline void mp_ins_vol_phase(const DG_FP *alpha, const DG_FP *s, DG_FP *vol) {
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    DG_FP step = tanh(PI * s[i] / *alpha);
    vol[i] = fmax(step, 0.0);
  }
}
