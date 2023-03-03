inline void p_multigrid_vec_normalise(const int *order, const DG_FP *in,
                                      const DG_FP *norm, DG_FP *out) {
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];
  for(int i = 0; i < dg_np; i++) {
    out[i] = in[i] / *norm;
  }
}
