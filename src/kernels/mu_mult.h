inline void mu_mult(const DG_FP *mu, DG_FP *out) {
  for(int i = 0; i < DG_NP; i++) {
    out[i] *= mu[i];
  }
}
