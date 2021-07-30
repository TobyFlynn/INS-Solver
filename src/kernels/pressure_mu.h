inline void pressure_mu(const double *mu, double *curl) {
  for(int i = 0; i < DG_NP; i++) {
    curl[i] *= mu[i];
  }
}
