inline void advec_2d_0(const double *q, const double *u, const double *v,
                    double *f, double *g) {
  for(int i = 0; i < DG_NP; i++) {
    f[i] = u[i] * q[i];
    g[i] = v[i] * q[i];
  }
}
