inline void grad_3d_geof_sp(const int *p, const DG_FP *geof, const float *ur,
                         const float *us, const float *ut, float *ux, float *uy,
                         float *uz) {
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
  for(int m = 0; m < dg_np; m++) {
    ux[m] = rx * ur[m] + sx * us[m] + tx * ut[m];
    uy[m] = ry * ur[m] + sy * us[m] + ty * ut[m];
    uz[m] = rz * ur[m] + sz * us[m] + tz * ut[m];
  }
}
