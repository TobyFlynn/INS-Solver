inline void p_multigrid_restriction(const int *p, const double *Au,
                                    const double *f, double *b, double *u) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * 5];
  for(int i = 0; i < dg_np; i++) {
    b[i] = f[i] - Au[i];
    u[i] = 0.0;
  }
}