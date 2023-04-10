inline void copy_diag(const int *order, const DG_FP *op, DG_FP *u) {
  // Au = u + (F - Au)
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];

  for(int i = 0; i < dg_np; i++) {
    const int op_ind = i * dg_np + i;
    u[i] = op[op_ind];
  }
}
