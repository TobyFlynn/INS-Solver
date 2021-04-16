inline void poisson_rhs_flux(const double *nx, const double *ny, const double *sJ,
                             const double *numFlux, double *fluxX, double *fluxY) {
  for(int i = 0; i < 21; i++) {
    fluxX[i] = nx[i] * gaussW_g[i % 7] * sJ[i] * numFlux[i];
    fluxY[i] = ny[i] * gaussW_g[i % 7] * sJ[i] * numFlux[i];
  }
}
