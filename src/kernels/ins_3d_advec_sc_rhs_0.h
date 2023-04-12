inline void ins_3d_advec_sc_rhs_0(const DG_FP *us, const DG_FP *vs, const DG_FP *ws,
                           const DG_FP *ub, const DG_FP *vb, const DG_FP *wb,
                           DG_FP *f00, DG_FP *f01, DG_FP *f02, DG_FP *f10,
                           DG_FP *f11, DG_FP *f12, DG_FP *f20, DG_FP *f21,
                           DG_FP *f22) {
  for(int i = 0; i < DG_NP; i++) {
    f00[i] = ub[i] * us[i]; f01[i] = ub[i] * vs[i]; f02[i] = ub[i] * ws[i];
    f10[i] = vb[i] * us[i]; f11[i] = vb[i] * vs[i]; f12[i] = vb[i] * ws[i];
    f20[i] = wb[i] * us[i]; f21[i] = wb[i] * vs[i]; f22[i] = wb[i] * ws[i];
  }
}
