inline void p_multigrid_relaxation_chebyshev_2(const int *order, const DG_FP *b,
                                               const DG_FP *diag, DG_FP *res) {
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];
  for(int i = 0; i < dg_np; i++) {
    res[i] -= b[i] / diag[i];
  }
}
