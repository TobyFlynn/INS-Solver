inline void pmf_3d_mult_cells_geof(const int *p, const DG_FP *geof, DG_FP *in_x,
                                   DG_FP *in_y, DG_FP *in_z) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  for(int n = 0; n < dg_np; n++) {
    const DG_FP x = in_x[n];
    const DG_FP y = in_y[n];
    const DG_FP z = in_z[n];
    in_x[n] = geof[RX_IND] * x + geof[RY_IND] * y + geof[RZ_IND] * z;
    in_y[n] = geof[SX_IND] * x + geof[SY_IND] * y + geof[SZ_IND] * z;
    in_z[n] = geof[TX_IND] * x + geof[TY_IND] * y + geof[TZ_IND] * z;
  }
}
