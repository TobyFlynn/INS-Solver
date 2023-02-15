inline void block_jacobi_pre(const int *p, const DG_FP *in, const DG_FP *pre, DG_FP *out) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  for(int i = 0; i < dg_np; i++) {
    out[i] = 0.0;
    for(int j = 0; j < dg_np; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      out[i] += pre[ind] * in[j];
    }
  }
}
