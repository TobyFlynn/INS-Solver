inline void ins_3d_advec_sc_rk_2(const DG_FP *dt, const DG_FP *u0, const DG_FP *v0,
                                 const DG_FP *w0, const DG_FP *u1, const DG_FP *v1,
                                 const DG_FP *w1, const DG_FP *u2, const DG_FP *v2,
                                 const DG_FP *w2, DG_FP *uT, DG_FP *vT, DG_FP *wT) {
  for(int i = 0; i < DG_NP; i++) {
    uT[i] -= *dt * (u0[i] / 6.0 + u1[i] / 6.0 + 2.0 * u2[i] / 3.0);
    vT[i] -= *dt * (v0[i] / 6.0 + v1[i] / 6.0 + 2.0 * v2[i] / 3.0);
    wT[i] -= *dt * (w0[i] / 6.0 + w1[i] / 6.0 + 2.0 * w2[i] / 3.0);
  }
}
