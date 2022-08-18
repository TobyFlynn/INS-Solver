inline void copy_dg_np(const double *in, double *out) {
  for(int i = 0; i < DG_NP; i++) {
    out[i] = in[i];
  }
}