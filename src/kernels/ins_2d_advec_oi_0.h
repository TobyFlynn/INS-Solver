inline void ins_2d_advec_oi_0(const DG_FP *geof, DG_FP *f00, DG_FP *f01, 
                              DG_FP *f10, DG_FP *f11) {
  for(int i = 0; i < DG_CUB_2D_NP; i++) {
    const DG_FP u = f00[i]; const DG_FP v = f01[i];
    // const DG_FP u = u_[i]; const DG_FP v = v_[i];
    const DG_FP tmp_dr = geof[RX_IND] * u + geof[RY_IND] * v;
    const DG_FP tmp_ds = geof[SX_IND] * u + geof[SY_IND] * v;
    f00[i] = u * tmp_dr; f01[i] = u * tmp_ds; 
    f10[i] = v * tmp_dr; f11[i] = v * tmp_ds; 
  }
}
