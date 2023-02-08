inline void advec_3d_0(const double *q, const double *u, const double *v,
                       const double *w, double *f, double *g, double *h) {
  for(int i = 0; i < DG_NP; i++) {
    f[i] = u[i] * q[i];
    g[i] = v[i] * q[i];
    h[i] = w[i] * q[i];
  }
}
