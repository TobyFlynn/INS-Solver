inline void ls_bcell_reset(const int *p, const double *x, const double *y, double *s) {
  // Get constants
  const int dg_np = DG_CONSTANTS[(*p - 1) * 5];

  for(int i = 0; i < dg_np; i++) {
      s[i] = sqrt((x[i] - 1.0) * (x[i] - 1.0) + (y[i] - 0.5) * (y[i] - 0.5)) - 0.2;
  }
}
