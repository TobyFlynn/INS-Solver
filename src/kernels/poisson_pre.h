inline void poisson_pre(const int *p, const double *in, const double *pre, double *out) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * 5];
  for(int i = 0; i < dg_np; i++) {
    out[i] = 0.0;
    for(int j = 0; j < dg_np; j++) {
      int ind = i + j * dg_np;
      out[i] += pre[ind] * in[j];
    }
  }
}
