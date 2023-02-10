inline void poisson_mm(const double *factor, const int *p, const double *mm,
                       double *op1) {
  // Get constants for this element's order
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  for(int m = 0; m < dg_np; m++) {
    for(int n = 0; n < dg_np; n++) {
      int ind = m + n * dg_np;
      op1[ind] += *factor * mm[ind];
    }
  }
}
