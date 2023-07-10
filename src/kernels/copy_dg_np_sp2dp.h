inline void copy_dg_np_sp2dp(const float *in, DG_FP *out) {
  for(int i = 0; i < DG_NP; i++) {
    out[i] = (double)in[i];
  }
}
