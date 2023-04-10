inline void p_multigrid_relaxation_chebyshev_0(const DG_FP *factor,
                                const int *order, const DG_FP *b,
                                const DG_FP *diag, DG_FP *res, DG_FP *out) {
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];
  for(int i = 0; i < dg_np; i++) {
    res[i] = b[i] - res[i];
    res[i] /= diag[i];
    out[i] = *factor * res[i];
  }
}
