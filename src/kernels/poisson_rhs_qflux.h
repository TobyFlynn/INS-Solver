inline void poisson_rhs_qflux(const double *nx, const double *ny, const double *sJ,
                              const double *tau, const double *gu, double *uFlux,
                              double *qxFlux, double *qyFlux, double *flux) {
  for(int i = 0; i < 21; i++) {
    // flux[i] = nx[i] * qxFlux[i] + ny[i] * qyFlux[i] - 2.0 * tau[i / 7] * (gu[i] - uFlux[i]);
    flux[i] = nx[i] * qxFlux[i] + ny[i] * qyFlux[i] - tau[i / 7] * (gu[i] - uFlux[i]);
    flux[i] *= gaussW_g[i % 7] * sJ[i];
  }

  for(int i = 0; i < 21; i++) {
    uFlux[i] = 0.0;
    qxFlux[i] = 0.0;
    qyFlux[i] = 0.0;
  }
}
