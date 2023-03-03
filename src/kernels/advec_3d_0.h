inline void advec_3d_0(const DG_FP *q, const DG_FP *u, const DG_FP *v,
                       const DG_FP *w, DG_FP *f, DG_FP *g, DG_FP *h) {
  for(int i = 0; i < DG_NP; i++) {
    f[i] = u[i] * q[i];
    g[i] = v[i] * q[i];
    h[i] = w[i] * q[i];
  }
}
