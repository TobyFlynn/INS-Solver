inline void poisson_mult_cells(const int *p, const DG_FP *u, const DG_FP *op,
                               DG_FP *rhs) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  for(int m = 0; m < dg_np; m++) {
    rhs[m] = 0.0;
    for(int n = 0; n < dg_np; n++) {
      int ind = m + n * dg_np;
      rhs[m] += op[ind] * u[n];
    }
  }
}
