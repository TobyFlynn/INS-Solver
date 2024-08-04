inline void fvmf_3d_grad(const DG_FP *geof, const DG_FP * fact, DG_FP *ux, DG_FP *uy,
                         DG_FP * uz, DG_FP *vx, DG_FP *vy, DG_FP *vz, DG_FP *wx, 
                         DG_FP *wy, DG_FP *wz) {
  const DG_FP rx = geof[RX_IND];
  const DG_FP sx = geof[SX_IND];
  const DG_FP tx = geof[TX_IND];
  const DG_FP ry = geof[RY_IND];
  const DG_FP sy = geof[SY_IND];
  const DG_FP ty = geof[TY_IND];
  const DG_FP rz = geof[RZ_IND];
  const DG_FP sz = geof[SZ_IND];
  const DG_FP tz = geof[TZ_IND];
  for(int m = 0; m < DG_NP; m++) {
    const DG_FP ur = ux[m];
    const DG_FP us = uy[m];
    const DG_FP ut = uz[m];
    const DG_FP _fact = fact[m];
    ux[m] = _fact * (rx * ur + sx * us + tx * ut);
    uy[m] = _fact * (ry * ur + sy * us + ty * ut);
    uz[m] = _fact * (rz * ur + sz * us + tz * ut);
    const DG_FP vr = vx[m];
    const DG_FP vs = vy[m];
    const DG_FP vt = vz[m];
    vx[m] = _fact * (rx * vr + sx * vs + tx * vt);
    vy[m] = _fact * (ry * vr + sy * vs + ty * vt);
    vz[m] = _fact * (rz * vr + sz * vs + tz * vt);
    const DG_FP wr = wx[m];
    const DG_FP ws = wy[m];
    const DG_FP wt = wz[m];
    wx[m] = _fact * (rx * wr + sx * ws + tx * wt);
    wy[m] = _fact * (ry * wr + sy * ws + ty * wt);
    wz[m] = _fact * (rz * wr + sz * ws + tz * wt);
  }
}