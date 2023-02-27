inline void block_jacobi_pre(const int *p, const DG_FP *in, const DG_FP *pre, DG_FP *out) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  op2_in_kernel_gemv(false, dg_np, dg_np, 1.0, pre, dg_np, in, 0.0, out);
}
