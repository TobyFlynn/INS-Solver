inline void vec_normalise(const int *p, const double *in, const double *norm, double *out) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  for(int i = 0; i < dg_np; i++) {
    out[i] = in[i] / *norm;
  }
}
