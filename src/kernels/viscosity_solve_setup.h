inline void viscosity_solve_setup(const double *mu, const double *rho,
                                  const double *mmConst, double *factor,
                                  double *mmFactor) {
  for(int i = 0; i < DG_NP; i++) {
    factor[i] = mu[i];
    mmFactor[i] = *mmConst * rho[i];
  }
}
