inline void ins_2d_advec_sc_rhs_oi_0(const DG_FP *geof, DG_FP *f00, DG_FP *f01,
                                     DG_FP *f10, DG_FP *f11) {
  for(int i = 0; i < DG_CUB_2D_NP; i++) {
    const DG_FP us = f00[i];
    const DG_FP vs = f01[i];
    const DG_FP ub = f10[i];
    const DG_FP vb = f11[i];
    const DG_FP tmp_dr = geof[RX_IND] * us + geof[RY_IND] * vs;
    const DG_FP tmp_ds = geof[SX_IND] * us + geof[SY_IND] * vs;
    f00[i] = ub * tmp_dr; f01[i] = ub * tmp_ds;
    f10[i] = vb * tmp_dr; f11[i] = vb * tmp_ds;
  }
}
