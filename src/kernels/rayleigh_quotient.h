inline void rayleigh_quotient(const int *p, const double *b, const double *Ab, double *top, double *bottom) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * 5];
  for(int i = 0; i < dg_np; i++) {
    *top += b[i] * Ab[i];
    *bottom += b[i] * b[i];
  }
}
