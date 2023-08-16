inline void ins_3d_advec_3(const DG_FP *a0, const DG_FP *a1, const DG_FP *b0,
                           const DG_FP *b1, const DG_FP *dt, const DG_FP *u,
                           const DG_FP *v, const DG_FP *w, const DG_FP *u_old,
                           const DG_FP *v_old, const DG_FP *w_old,
                           const DG_FP *n0, const DG_FP *n1, const DG_FP *n2,
                           const DG_FP *n0_old, const DG_FP *n1_old,
                           const DG_FP *n2_old, DG_FP *uT, DG_FP *vT,
                           DG_FP *wT) {
  const DG_FP _a0 = *a0;
  const DG_FP _a1 = *a1;
  const DG_FP _b0 = *b0;
  const DG_FP _b1 = *b1;
  const DG_FP _dt = *dt;
  for(int i = 0; i < DG_NP; i++) {
    uT[i] = (_a0 * u[i] + _a1 * u_old[i]) - _dt * (_b0 * n0[i] + _b1 * n0_old[i]);
    vT[i] = (_a0 * v[i] + _a1 * v_old[i]) - _dt * (_b0 * n1[i] + _b1 * n1_old[i]);
    wT[i] = (_a0 * w[i] + _a1 * w_old[i]) - _dt * (_b0 * n2[i] + _b1 * n2_old[i]);
  }
}
