inline void mp_ins_sv_phase(const DG_FP *alpha, const DG_FP *s,
                            const DG_FP *v, DG_FP *sv) {
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    DG_FP step = tanh(PI * s[i] / *alpha);
    sv[i] = fmax(step, 0.0) * v[i];
  }
}
