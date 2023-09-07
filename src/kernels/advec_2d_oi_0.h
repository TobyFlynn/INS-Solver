inline void advec_2d_oi_0(const DG_FP *geof, const DG_FP *q, DG_FP *f, DG_FP *g) {
  for(int i = 0; i < DG_NP; i++) {
    const DG_FP u = f[i];
    const DG_FP v = g[i];
    const DG_FP _f = u * q[i];
    const DG_FP _g = v * q[i];
    f[i] = geof[RX_IND] * _f + geof[RY_IND] * _g;
    g[i] = geof[SX_IND] * _f + geof[SY_IND] * _g;
  }
}
