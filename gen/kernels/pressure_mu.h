inline void pressure_mu(const double *mu, double *curl) {
  for(int i = 0; i < 6; i++) {
    curl[i] *= mu[i];
  }
}
