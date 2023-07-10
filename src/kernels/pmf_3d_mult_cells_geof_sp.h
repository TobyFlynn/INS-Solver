inline void pmf_3d_mult_cells_geof_sp(const int *p, const DG_FP *geof, float *in_x,
                                   float *in_y, float *in_z) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  const float rx = (float)geof[RX_IND];
  const float sx = (float)geof[SX_IND];
  const float tx = (float)geof[TX_IND];
  const float ry = (float)geof[RY_IND];
  const float sy = (float)geof[SY_IND];
  const float ty = (float)geof[TY_IND];
  const float rz = (float)geof[RZ_IND];
  const float sz = (float)geof[SZ_IND];
  const float tz = (float)geof[TZ_IND];
  for(int n = 0; n < dg_np; n++) {
    const float x = in_x[n];
    const float y = in_y[n];
    const float z = in_z[n];
    in_x[n] = rx * x + ry * y + rz * z;
    in_y[n] = sx * x + sy * y + sz * z;
    in_z[n] = tx * x + ty * y + tz * z;
  }
}
