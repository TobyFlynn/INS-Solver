inline void ins_3d_advec_sc_copy(const DG_FP *u0, const DG_FP *v0, const DG_FP *w0,
                                 DG_FP *u1, DG_FP *v1, DG_FP *w1) {
  for(int i = 0; i < DG_NP; i++) {
    u1[i] = u0[i];
    v1[i] = v0[i];
    w1[i] = w0[i];
  }
}
