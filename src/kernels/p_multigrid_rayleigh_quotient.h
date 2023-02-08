inline void p_multigrid_rayleigh_quotient(const int *order, const double *b,
                                          const double *Ab, double *top,
                                          double *bottom) {
  const int dg_np = DG_CONSTANTS[(*order - 1) * 5];
  for(int i = 0; i < dg_np; i++) {
    *top += b[i] * Ab[i];
    *bottom += b[i] * b[i];
  }
}
