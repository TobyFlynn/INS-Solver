inline void poisson_set_rhs1(const double *nx, const double *ny, const double *sJ,
                             const double *bc, double *bcNX, double *bcNY, double *tau) {
  for(int i = 0; i < 15; i++) {
    bcNX[i] = sJ[i] * nx[i] * bc[i];
    bcNY[i] = sJ[i] * ny[i] * bc[i];
    tau[i]  = sJ[i] * tau[i];
  }
}
