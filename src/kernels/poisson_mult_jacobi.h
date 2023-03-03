inline void poisson_mult_jacobi(const int *p, const DG_FP *op, DG_FP *rhs) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  for(int i = 0; i < dg_np; i++) {
    // const int op_ind = i * dg_np + i;
    int op_ind = DG_MAT_IND(i, i, dg_np, dg_np);
    rhs[i] = rhs[i] / op[op_ind];
  }
}
