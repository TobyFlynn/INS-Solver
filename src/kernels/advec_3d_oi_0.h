inline void advec_3d_oi_0(const DG_FP *geof, const DG_FP *q, DG_FP *f, DG_FP *g, DG_FP *h) {
  for(int i = 0; i < DG_CUB_3D_NP; i++) {
    const DG_FP u = f[i];
    const DG_FP v = g[i];
    const DG_FP w = h[i];
    const DG_FP _f = u * q[i];
    const DG_FP _g = v * q[i];
    const DG_FP _h = w * q[i];
    f[i] = geof[RX_IND] * _f + geof[RY_IND] * _g + geof[RZ_IND] * _h;
    g[i] = geof[SX_IND] * _f + geof[SY_IND] * _g + geof[SZ_IND] * _h;
    h[i] = geof[TX_IND] * _f + geof[TY_IND] * _g + geof[TZ_IND] * _h;
  }
}
