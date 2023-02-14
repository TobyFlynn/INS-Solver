inline void advec_2d_0(const DG_FP *q, const DG_FP *u, const DG_FP *v,
                    DG_FP *f, DG_FP *g) {
  for(int i = 0; i < DG_NP; i++) {
    f[i] = u[i] * q[i];
    g[i] = v[i] * q[i];
  }
}
