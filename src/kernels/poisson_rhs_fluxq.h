inline void poisson_rhs_fluxq(const double *nx, const double *ny, const double *fscale,
                       const double *tau, const double *u, const double *du, const double *qx,
                       const double *qy, double *exQx, double *exQy, double *fluxq) {
  for(int i = 0; i < 15; i++) {
    double dqx = (qx[FMASK[i]] + exQx[i]) / 2.0;
    double dqy = (qy[FMASK[i]] + exQy[i]) / 2.0;
    fluxq[i] = fscale[i] * (nx[i] * dqx + ny[i] * dqy + tau[i] * (u[FMASK[i]] - du[i]));
    exQx[i] = 0.0;
    exQy[i] = 0.0;
  }
}
