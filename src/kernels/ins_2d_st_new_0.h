inline void ins_2d_st_new_0(const DG_FP *alpha_, const DG_FP *s, DG_FP *out) {
  const DG_FP alpha = *alpha_;
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    DG_FP step = 0.5 * tanh(PI * s[i] / alpha);
    out[i] = step;
  }
}