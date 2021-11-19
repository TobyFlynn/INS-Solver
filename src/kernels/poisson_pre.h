inline void poisson_pre(const double *in, const double *pre, double *out) {
  for(int i = 0; i < DG_NP; i++) {
    out[i] = 0.0;
    for(int j = 0; j < DG_NP; j++) {
      int ind = i + j * DG_NP;
      out[i] += pre[ind] * in[j];
    }
  }
}
