inline void poisson_pre(const double *in, const double *pre, double *out) {
  for(int i = 0; i < 6; i++) {
    out[i] = 0.0;
    for(int j = 0; j < 6; j++) {
      // int mm_ind = j * 15 + i;
      int ind = i * 6 + j;
      out[i] += pre[ind] * in[j];
    }
  }
}
