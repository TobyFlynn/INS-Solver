inline void pRHS_fluxq(const double *nx, const double *ny, const double *tau,
                       const double *du, const double *qx, const double *qy,
                       double *exQx, double *exQy, double *fluxq) {
  for(int i = 0; i < 15; i++) {
    double dqx = qx[FMASK[i]] - exQx[i];
    double dqy = qy[FMASK[i]] - exQy[i];
    fluxq[i] = (nx[i] * dqx + ny[i] * dqy + tau[i] * du[i]) / 2.0;
    exQx[i] = 0.0;
    exQy[i] = 0.0;
  }
}
