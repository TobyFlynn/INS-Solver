inline void ls_step_mu_vol_oi(const DG_FP *alpha, DG_FP *mu) {
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_CUB_2D_NP; i++) {
    DG_FP step = tanh(PI * mu[i] / *alpha);
    mu[i] = 0.5 * mu0 * (1.0 + step) + 0.5 * mu1 * (1.0 - step);
  }
}
