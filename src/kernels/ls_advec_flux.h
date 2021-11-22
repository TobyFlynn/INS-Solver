inline void ls_advec_flux(const double *q, const double *u, const double *v,
                          double *F, double *G) {
  for(int i = 0; i < DG_NP; i++) {
    F[i] = u[i] * q[i];
    G[i] = v[i] * q[i];
  }
}
