inline void ls_advec_rhs(const int *p, const double *dFdr, const double *dFds,
                         const double *dGdr, const double *dGds,
                         const double *rx, const double *ry, const double *sx,
                         const double *sy, const double *q, const double *exQ,
                         const double *u, const double *v, const double *sJ,
                         const double *nx, const double *ny, double *nFlux,
                         double *output) {
  // Get constants
  // const int dg_np  = DG_CONSTANTS[(*p - 1) * 5];
  const int dg_npf = DG_CONSTANTS[(*p - 1) * 5 + 1];
  const double *gaussW = &gaussW_g[(*p - 1) * DG_GF_NP];

  for(int i = 0; i < DG_NP; i++) {
    // output[i] = rx[i] * dFdr[i] + sx[i] * dFds[i] + ry[i] * dGdr[i] + sy[i] * dGds[i];
  }

  double mQ[DG_G_NP];
  double mF[DG_G_NP];
  double mG[DG_G_NP];
  for(int i = 0; i < DG_G_NP; i++) {
    mQ[i] = q[i];
    mF[i] = u[i] * q[i];
    mG[i] = v[i] * q[i];
  }

  double pF[DG_G_NP];
  double pG[DG_G_NP];
  for(int i = 0; i < DG_G_NP; i++) {
    pF[i]  = u[i] * exQ[i];
    pG[i]  = v[i] * exQ[i];
  }

  // Lift fluxes
  for(int i = 0; i < DG_G_NP; i++) {
    // double a = sqrt(u[ind] * u[ind] + v[ind] * v[ind]);
    // nFlux[i] = nx[i] * (mF[i] + pF[i]) + ny[i] * (mG[i] + pG[i]) + a * (mQ[i] - exQ[i]);
    // nFlux[i] *= 0.5 * fscale[i];

    // nFlux[i] = (nx[i] * u[i] + ny[i] * v[i]) * (q[i] + exQ[i]);
    // nFlux[i] += fabs(nx[i] * u[i] + ny[i] * v[i]) * (q[i] - exQ[i]);
    // nFlux[i] *= 0.5 * sJ[i] * gaussW[i % DG_GF_NP];

    double a = sqrt(u[i] * u[i] + v[i] * v[i]);
    nFlux[i] = nx[i] * (mF[i] + pF[i]) + ny[i] * (mG[i] + pG[i]) + a * (mQ[i] - exQ[i]);
    nFlux[i] *= 0.5 * sJ[i] * gaussW[i % DG_GF_NP];
  }
}
