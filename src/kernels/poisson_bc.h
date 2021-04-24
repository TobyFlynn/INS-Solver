inline void poisson_bc(const double *nx, const double *ny, const double *sJ,
                       const double *bc, double *gx, double *gy) {
  for(int i = 0; i < 21; i++) {
    gx[i] = gaussW_g[i % 7] * sJ[i] * nx[i] * bc[i];
    gy[i] = gaussW_g[i % 7] * sJ[i] * ny[i] * bc[i];
  }
}
