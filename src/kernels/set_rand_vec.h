inline void set_rand_vec(const DG_FP *in, DG_FP *out) {
  for(int i = 0; i < DG_NP; i++) {
    out[i] = in[i];
  }
}