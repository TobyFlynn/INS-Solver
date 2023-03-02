inline void poisson_mult_cells_coarse(const DG_FP *u, const DG_FP *op, DG_FP *rhs) {
  op2_in_kernel_gemv(false, DG_NP_N1, DG_NP_N1, 1.0, op, DG_NP_N1, u, 0.0, rhs);
}
