inline void pmf_2d_mult_cells_geof(const int *p, const DG_FP *geof, DG_FP *in_x,
                                   DG_FP *in_y) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  const DG_FP rx = geof[RX_IND];
  const DG_FP sx = geof[SX_IND];
  const DG_FP ry = geof[RY_IND];
  const DG_FP sy = geof[SY_IND];
  for(int n = 0; n < dg_np; n++) {
    const DG_FP x = in_x[n];
    const DG_FP y = in_y[n];
    in_x[n] = rx * x + ry * y;
    in_y[n] = sx * x + sy * y;
  }
}
