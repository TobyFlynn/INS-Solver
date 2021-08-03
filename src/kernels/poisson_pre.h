inline void poisson_pre(const double *in, const double *pre, double *out) {
  for(int i = 0; i < DG_NP; i++) {
    out[i] = 0.0;
    for(int j = 0; j < DG_NP; j++) {
      // int mm_ind = j * 15 + i;
      int ind = i * DG_NP + j;
      out[i] += pre[ind] * in[j];
    }
  }
}
