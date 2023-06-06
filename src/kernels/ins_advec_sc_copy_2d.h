inline void ins_advec_sc_copy_2d(const DG_FP *u0, const DG_FP *v0, DG_FP *u1,
                                 DG_FP *v1) {
  for(int i = 0; i < DG_NP; i++) {
    u1[i] = u0[i];
    v1[i] = v0[i];
  }
}
