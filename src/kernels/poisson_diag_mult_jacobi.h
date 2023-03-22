inline void poisson_diag_mult_jacobi(const int *p, const DG_FP *diag, DG_FP *rhs) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  for(int i = 0; i < dg_np; i++) {
    rhs[i] = rhs[i] / diag[i];
  }
}
