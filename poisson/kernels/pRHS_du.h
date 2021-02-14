inline void pRHS_du(const double *nx, const double *ny, const double *fscale,
                    const double *U, double *exU, double *du,
                    double *fluxXu, double *fluxYu) {
  for(int i = 0; i < 15; i++) {
    du[i] = U[i] - exU[i];
    fluxXu[i] = fscale[i] * (nx[i] * du[i] / 2.0);
    fluxYu[i] = fscale[i] * (ny[i] * du[i] / 2.0);
    exU[i] = 0.0;
  }
}
