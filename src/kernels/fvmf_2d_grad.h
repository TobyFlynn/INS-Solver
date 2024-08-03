inline void fvmf_2d_grad(const DG_FP *geof, const DG_FP *factor,
                         DG_FP *ux, DG_FP *uy,DG_FP *vx, DG_FP *vy) {
  for(int m = 0; m < DG_NP; m++) {
    const DG_FP ur = ux[m];
    const DG_FP us = uy[m];
    const DG_FP vr = vx[m];
    const DG_FP vs = vy[m];
    ux[m] = factor[m] * (geof[RX_IND] * ur + geof[SX_IND] * us);
    uy[m] = factor[m] * (geof[RY_IND] * ur + geof[SY_IND] * us);
    vx[m] = factor[m] * (geof[RX_IND] * vr + geof[SX_IND] * vs);
    vy[m] = factor[m] * (geof[RY_IND] * vr + geof[SY_IND] * vs);
  }
}