inline void ls_advec_rhs(const int *p, const double *dFdr, const double *dFds,
                         const double *dGdr, const double *dGds,
                         const double *rx, const double *ry, const double *sx,
                         const double *sy, const double *q, const double *exQ,
                         const double *u, const double *v, const double *fscale,
                         const double *nx, const double *ny, double *nFlux,
                         double *output) {
  // Get constants
  const int dg_np  = DG_CONSTANTS[(*p - 1) * 5];
  const int dg_npf = DG_CONSTANTS[(*p - 1) * 5 + 1];
  const int *fmask = &FMASK[(*p - 1) * 3 * DG_NPF];

  for(int i = 0; i < dg_np; i++) {
    output[i] = rx[i] * dFdr[i] + sx[i] * dFds[i] + ry[i] * dGdr[i] + sy[i] * dGds[i];
  }

  double mQ[3 * DG_NPF];
  double mF[3 * DG_NPF];
  double mG[3 * DG_NPF];
  for(int i = 0; i < 3 * dg_npf; i++) {
    int ind = fmask[i];
    mQ[i] = q[ind];
    mF[i] = u[ind] * q[ind];
    mG[i] = v[ind] * q[ind];
  }

  double pF[3 * DG_NPF];
  double pG[3 * DG_NPF];
  for(int i = 0; i < 3 * dg_npf; i++) {
    int ind = fmask[i];
    pF[i]  = u[ind] * exQ[i];
    pG[i]  = v[ind] * exQ[i];
  }

  // Lift fluxes
  for(int i = 0; i < 3 * dg_npf; i++) {
    int ind = fmask[i];
    // double a = sqrt(u[ind] * u[ind] + v[ind] * v[ind]);
    // nFlux[i] = nx[i] * (mF[i] + pF[i]) + ny[i] * (mG[i] + pG[i]) + a * (mQ[i] - exQ[i]);
    // nFlux[i] *= 0.5 * fscale[i];
    nFlux[i] = (nx[i] * u[ind] + ny[i] * v[ind]) * (q[ind] + exQ[i]);
    nFlux[i] += fabs(nx[i] * u[ind] + ny[i] * v[ind]) * (q[ind] - exQ[i]);
    nFlux[i] *= 0.5 * fscale[i];
  }
}
