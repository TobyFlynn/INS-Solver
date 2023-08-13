inline void ins_3d_advec_oi_0(const DG_FP *geof, DG_FP *f00, DG_FP *f01, 
                    DG_FP *f02, DG_FP *f10, DG_FP *f11, DG_FP *f12, 
                    DG_FP *f20, DG_FP *f21, DG_FP *f22) {
  for(int i = 0; i < DG_CUB_3D_NP; i++) {
    const DG_FP u = f00[i]; const DG_FP v = f01[i]; const DG_FP w = f02[i];
    const DG_FP tmp_dr = geof[RX_IND] * u + geof[RY_IND] * v + geof[RZ_IND] * w;
    const DG_FP tmp_ds = geof[SX_IND] * u + geof[SY_IND] * v + geof[SZ_IND] * w;
    const DG_FP tmp_dt = geof[TX_IND] * u + geof[TY_IND] * v + geof[TZ_IND] * w;
    f00[i] = u * tmp_dr; f01[i] = u * tmp_ds; f02[i] = u * tmp_dt;
    f10[i] = v * tmp_dr; f11[i] = v * tmp_ds; f12[i] = v * tmp_dt;
    f20[i] = w * tmp_dr; f21[i] = w * tmp_ds; f22[i] = w * tmp_dt;
  }
}
