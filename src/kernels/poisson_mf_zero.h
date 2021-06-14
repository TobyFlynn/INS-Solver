inline void poisson_mf_zero(double *flux0, double *flux1, double *flux2) {
  for(int i = 0; i < 21; i++) {
    flux0[i] = 0.0;
    flux1[i] = 0.0;
    flux2[i] = 0.0;
  }
}
