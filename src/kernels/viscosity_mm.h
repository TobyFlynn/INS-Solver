inline void viscosity_mm(const int *p, const double *mm, const double *v,
                         double *res) {
  // Get constants
  const int dg_np = DG_CONSTANTS[(*p - 1) * 5];
  for(int m = 0; m < dg_np; m++) {
    res[m] = 0.0;
    for(int n = 0; n < dg_np; n++) {
      int ind = m + n * dg_np;
      res[m] += mm[ind] * v[n];
    }
  }
}
