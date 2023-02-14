inline void poisson_vis_fact(const DG_FP *mmConst, const DG_FP *mu,
                             const DG_FP *rho, DG_FP *factor,
                             DG_FP *mmFactor) {
  for(int i = 0; i < DG_NP; i++) {
    factor[i]   = mu[i];
    mmFactor[i] = *mmConst * rho[i];
  }
}
