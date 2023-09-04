inline void ins_2d_st_1(const DG_FP *geof, DG_FP *dx, DG_FP *dy) {
  for(int i = 0; i < DG_NP; i++) {
    const DG_FP dr = dx[i];
    const DG_FP ds = dy[i];
    dx[i] = geof[RX_IND] * dr + geof[SX_IND] * ds;
    dy[i] = geof[RY_IND] * dr + geof[SY_IND] * ds;
  }
}