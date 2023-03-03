inline void p_multigrid_rayleigh_quotient(const int *order, const DG_FP *b,
                                          const DG_FP *Ab, DG_FP *top,
                                          DG_FP *bottom) {
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];
  for(int i = 0; i < dg_np; i++) {
    *top += b[i] * Ab[i];
    *bottom += b[i] * b[i];
  }
}
