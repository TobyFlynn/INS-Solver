inline void ins_3d_st_1(const DG_FP *geof, DG_FP *dx, DG_FP *dy, DG_FP *dz) {
  const DG_FP rx = geof[RX_IND];
  const DG_FP sx = geof[SX_IND];
  const DG_FP tx = geof[TX_IND];
  const DG_FP ry = geof[RY_IND];
  const DG_FP sy = geof[SY_IND];
  const DG_FP ty = geof[TY_IND];
  const DG_FP rz = geof[RZ_IND];
  const DG_FP sz = geof[SZ_IND];
  const DG_FP tz = geof[TZ_IND];
  for(int i = 0; i < DG_NP; i++) {
    const DG_FP dr = dx[i];
    const DG_FP ds = dy[i];
    const DG_FP dt = dz[i];
    dx[i] = rx * dr + sx * ds + tx * dt;
    dy[i] = ry * dr + sy * ds + ty * dt;
    dz[i] = rz * dr + sz * ds + tz * dt;
  }
}
