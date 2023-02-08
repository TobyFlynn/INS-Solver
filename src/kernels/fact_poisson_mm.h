inline void fact_poisson_mm(const int *p, const double *mm, const double *factor,
                      double *op1) {
  // Get constants for this element's order
  const int dg_np = DG_CONSTANTS[(*p - 1) * 5];

  for(int m = 0; m < dg_np; m++) {
    for(int n = 0; n < dg_np; n++) {
      int ind = m + n * dg_np;
      op1[ind] += factor[m] * mm[ind];
    }
  }
}
