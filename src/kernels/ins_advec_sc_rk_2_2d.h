inline void ins_advec_sc_rk_2_2d(const DG_FP *dt, const DG_FP *u0, const DG_FP *v0,
                                 const DG_FP *u1, const DG_FP *v1, DG_FP *uT,
                                 DG_FP *vT) {
  for(int i = 0; i < DG_NP; i++) {
    uT[i] -= *dt * (u0[i] + 2.0 * u1[i] / 3.0);
    vT[i] -= *dt * (v0[i] + 2.0 * v1[i] / 3.0);
  }
}
