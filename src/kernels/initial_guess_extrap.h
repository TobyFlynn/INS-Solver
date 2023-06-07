inline void initial_guess_extrap(const DG_FP *beta, const DG_FP *h0,
                                 const DG_FP *h1, const DG_FP *h2,
                                 const DG_FP *h3, DG_FP *out) {
  for(int i = 0; i < DG_NP; i++) {
    out[i] = beta[0] * h0[i] + beta[1] * h1[i] + beta[2] * h2[i] + beta[3] * h3[i];
  }
}
