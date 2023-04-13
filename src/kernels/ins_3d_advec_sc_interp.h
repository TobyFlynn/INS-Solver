inline void ins_3d_advec_sc_interp(const DG_FP *t0, const DG_FP *t1, const DG_FP *tI,
                           const DG_FP *u0, const DG_FP *v0, const DG_FP *w0,
                           const DG_FP *u1, const DG_FP *v1, const DG_FP *w1,
                           DG_FP *uI, DG_FP *vI, DG_FP *wI) {
  const DG_FP diff_t = *t1 - *t0;
  const DG_FP dt = (*tI - *t1) / diff_t;
  for(int i = 0; i < DG_NP; i++) {
    const DG_FP diff_u = u1[i] - u0[i];
    const DG_FP diff_v = v1[i] - v0[i];
    const DG_FP diff_w = w1[i] - w0[i];
    uI[i] = u1[i] + diff_u * dt;
    vI[i] = v1[i] + diff_v * dt;
    wI[i] = w1[i] + diff_w * dt;
  }
}
