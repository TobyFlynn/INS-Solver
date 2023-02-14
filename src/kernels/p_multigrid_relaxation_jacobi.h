inline void p_multigrid_relaxation_jacobi(const DG_FP *dt, const int *order,
                                          const DG_FP *Au, const DG_FP *f,
                                          const DG_FP *op, DG_FP *u) {
  // Au = u + (F - Au)
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];

  for(int i = 0; i < dg_np; i++) {
    const int op_ind = i * dg_np + i;
    u[i] += *dt * (f[i] - Au[i]) / op[op_ind];
  }
}
