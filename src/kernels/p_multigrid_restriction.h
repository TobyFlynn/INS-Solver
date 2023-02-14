inline void p_multigrid_restriction(const int *p, const DG_FP *Au,
                                    const DG_FP *f, DG_FP *b, DG_FP *u) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  for(int i = 0; i < dg_np; i++) {
    b[i] = f[i] - Au[i];
    u[i] = 0.0;
  }
}
