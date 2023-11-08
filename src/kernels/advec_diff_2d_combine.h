inline void advec_diff_2d_combine(const DG_FP *tmp, DG_FP *out) {
  for(int i = 0; i < DG_NP; i++) {
    out[i] += tmp[i];
  }
}