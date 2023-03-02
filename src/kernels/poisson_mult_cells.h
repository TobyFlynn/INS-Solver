inline void poisson_mult_cells(const int *p, const DG_FP *u, const DG_FP *op,
                               DG_FP *rhs) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  op2_in_kernel_gemv(false, dg_np, dg_np, 1.0, op, dg_np, u, 0.0, rhs);
}
