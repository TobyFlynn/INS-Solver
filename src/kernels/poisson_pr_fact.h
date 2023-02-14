inline void poisson_pr_fact(const DG_FP *rho, DG_FP *factor) {
  for(int i = 0; i < DG_NP; i++) {
    factor[i] = 1.0 / rho[i];
  }
}
