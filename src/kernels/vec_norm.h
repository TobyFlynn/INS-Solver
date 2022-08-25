inline void vec_norm(const int *p, const double *vec, double *norm) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * 5];
  for(int i = 0; i < dg_np; i++) {
    *norm += vec[i] * vec[i];
  }
}