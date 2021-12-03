inline void poisson_vis_fact(const double *mmConst, const double *mu,
                             const double *rho, double *factor,
                             double *mmFactor) {
  for(int i = 0; i < DG_NP; i++) {
    factor[i]   = mu[i];
    mmFactor[i] = *mmConst * rho[i];
  }
}
