inline void ins_advec_sc_interp_2d(const DG_FP *t0, const DG_FP *t1,
                           const DG_FP *tI, const DG_FP *u0, const DG_FP *v0,
                           const DG_FP *u1, const DG_FP *v1, DG_FP *uI,
                           DG_FP *vI) {
  const DG_FP diff_t = *t1 - *t0;
  const DG_FP dt = (*tI - *t1) / diff_t;
  for(int i = 0; i < DG_NP; i++) {
    const DG_FP diff_u = u1[i] - u0[i];
    const DG_FP diff_v = v1[i] - v0[i];
    uI[i] = u1[i] + diff_u * dt;
    vI[i] = v1[i] + diff_v * dt;
  }
}
