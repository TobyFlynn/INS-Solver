inline void ins_2d_advec_sc_force_2d(const DG_FP *b0, const DG_FP *b1,
                              const DG_FP *dt, const DG_FP *F0, const DG_FP *F1,
                              const DG_FP *F0Old, const DG_FP *F1Old, DG_FP *uT,
                              DG_FP *vT) {
  for(int i = 0; i < DG_NP; i++) {
    uT[i] -= *dt * (*b0 * F0[i] + *b1 * F0Old[i]);
    vT[i] -= *dt * (*b0 * F1[i] + *b1 * F1Old[i]);
  }
}
