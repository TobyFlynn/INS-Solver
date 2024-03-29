inline void ins_3d_advec_sc_update(const DG_FP *a0, const DG_FP *a1,
                                   const DG_FP *u0, const DG_FP *v0, const DG_FP *w0,
                                   DG_FP *u1, DG_FP *v1, DG_FP *w1) {
  for(int i = 0; i < DG_NP; i++) {
    u1[i] = *a0 * u0[i] + *a1 * u1[i];
    v1[i] = *a0 * v0[i] + *a1 * v1[i];
    w1[i] = *a0 * w0[i] + *a1 * w1[i];
  }
}
