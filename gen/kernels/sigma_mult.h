inline void sigma_mult(const double *vis, const double *s, double *out) {
  double k = sqrt(*vis);
  for(int i = 0; i < 10; i++) {
    out[i] = k * s[i];
  }
}
