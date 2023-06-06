inline void ins_advec_sc_rhs_0_2d(const DG_FP *us, const DG_FP *vs,
                           const DG_FP *ub, const DG_FP *vb, DG_FP *f00,
                           DG_FP *f01, DG_FP *f10, DG_FP *f11) {
  for(int i = 0; i < DG_NP; i++) {
    f00[i] = ub[i] * us[i]; f01[i] = ub[i] * vs[i];
    f10[i] = vb[i] * us[i]; f11[i] = vb[i] * vs[i];
  }
}
