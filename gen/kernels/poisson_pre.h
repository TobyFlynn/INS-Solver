inline void poisson_pre(const double *in, const double *pre, double *out) {
  for(int i = 0; i < 10; i++) {
    out[i] = 0.0;
    for(int j = 0; j < 10; j++) {
      // int mm_ind = j * 15 + i;
      int ind = i * 10 + j;
      out[i] += pre[ind] * in[j];
    }
  }
}