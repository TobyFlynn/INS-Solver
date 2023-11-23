inline void reciprocal(const DG_FP *in, DG_FP *out) {
  for(int i = 0; i < DG_NP; i++) {
    out[i] = 1.0 / in[i];
  }
}