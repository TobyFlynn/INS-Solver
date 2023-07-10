inline void copy_dg_np_dp2sp(const DG_FP *in, float *out) {
  for(int i = 0; i < DG_NP; i++) {
    out[i] = (float)in[i];
  }
}
