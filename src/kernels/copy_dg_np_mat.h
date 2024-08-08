inline void copy_dg_np_mat(const DG_FP *in, DG_FP *out) {
  for(int i = 0; i < DG_NP * DG_NP; i++) {
    out[i] = in[i];
  }
}