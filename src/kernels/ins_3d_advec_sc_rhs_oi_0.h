inline void ins_3d_advec_sc_rhs_oi_0(const DG_FP *geof, DG_FP *f00, DG_FP *f01, 
                          DG_FP *f02, DG_FP *f10, DG_FP *f11, DG_FP *f12, 
                          DG_FP *f20, DG_FP *f21, DG_FP *f22) {
  for(int i = 0; i < DG_CUB_3D_NP; i++) {
    const DG_FP us = f00[i];
    const DG_FP vs = f01[i];
    const DG_FP ws = f02[i];
    const DG_FP ub = f10[i];
    const DG_FP vb = f11[i];
    const DG_FP wb = f12[i];
    const DG_FP tmp_dr = geof[RX_IND] * us + geof[RY_IND] * vs + geof[RZ_IND] * ws;
    const DG_FP tmp_ds = geof[SX_IND] * us + geof[SY_IND] * vs + geof[SZ_IND] * ws;
    const DG_FP tmp_dt = geof[TX_IND] * us + geof[TY_IND] * vs + geof[TZ_IND] * ws;
    f00[i] = ub * tmp_dr; f01[i] = ub * tmp_ds; f02[i] = ub * tmp_dt;
    f10[i] = vb * tmp_dr; f11[i] = vb * tmp_ds; f12[i] = vb * tmp_dt;
    f20[i] = wb * tmp_dr; f21[i] = wb * tmp_ds; f22[i] = wb * tmp_dt;
  }
}
