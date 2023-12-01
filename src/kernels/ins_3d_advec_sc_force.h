inline void ins_3d_advec_sc_force(const DG_FP *b0, const DG_FP *b1, const DG_FP *dt,
                              const DG_FP *f0, const DG_FP *f1, const DG_FP *f2,
                              const DG_FP *f0_old, const DG_FP *f1_old,
                              const DG_FP *f2_old, DG_FP *uT, DG_FP *vT,
                              DG_FP *wT) {
  const DG_FP _b0 = *b0;
  const DG_FP _b1 = *b1;
  const DG_FP _dt = *dt;
  for(int i = 0; i < DG_NP; i++) {
    uT[i] -= _dt * (_b0 * f0[i] + _b1 * f0_old[i]);
    vT[i] -= _dt * (_b0 * f1[i] + _b1 * f1_old[i]);
    wT[i] -= _dt * (_b0 * f2[i] + _b1 * f2_old[i]);
  }
}
