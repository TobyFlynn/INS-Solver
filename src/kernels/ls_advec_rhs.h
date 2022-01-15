inline void ls_advec_rhs(const int *p, const double *dFdr, const double *dFds,
                         const double *dGdr, const double *dGds,
                         const double *rx, const double *ry, const double *sx,
                         const double *sy, const double *q, const double *exQ,
                         const double *u, const double *v, const double *sJ,
                         const double *nx, const double *ny, double *nFlux,
                         double *output) {
  // Get constants
  // const int dg_np  = DG_CONSTANTS[(*p - 1) * 5];
  const double *gaussW = &gaussW_g[(*p - 1) * DG_GF_NP];

  for(int i = 0; i < DG_NP; i++) {
    output[i] = rx[i] * dFdr[i] + sx[i] * dFds[i] + ry[i] * dGdr[i] + sy[i] * dGds[i];
  }

  // Lift fluxes
  for(int i = 0; i < DG_G_NP; i++) {
    nFlux[i] = (nx[i] * u[i] + ny[i] * v[i]) * (q[i] + exQ[i]);
    nFlux[i] += fabs(nx[i] * u[i] + ny[i] * v[i]) * (q[i] - exQ[i]);
    nFlux[i] *= 0.5 * sJ[i] * gaussW[i % DG_GF_NP];
  }
}
