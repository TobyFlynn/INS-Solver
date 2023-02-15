inline void poisson_mm(const DG_FP *factor, const int *p, const DG_FP *mm,
                       DG_FP *op1) {
  // Get constants for this element's order
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  for(int m = 0; m < dg_np; m++) {
    for(int n = 0; n < dg_np; n++) {
      // int ind = m + n * dg_np;
      int ind = DG_MAT_IND(m, n, dg_np, dg_np);
      op1[ind] += *factor * mm[ind];
    }
  }
}
