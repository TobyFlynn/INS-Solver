inline void vmf_2d_mult_cells_geof(const DG_FP *geof, DG_FP *u_x, DG_FP *u_y, 
                                   DG_FP *v_x, DG_FP *v_y) {
  const DG_FP rx = geof[RX_IND];
  const DG_FP sx = geof[SX_IND];
  const DG_FP ry = geof[RY_IND];
  const DG_FP sy = geof[SY_IND];
  for(int n = 0; n < DG_NP; n++) {
    const DG_FP x = u_x[n];
    const DG_FP y = u_y[n];
    u_x[n] = rx * x + ry * y;
    u_y[n] = sx * x + sy * y;
  }

  for(int n = 0; n < DG_NP; n++) {
    const DG_FP x = v_x[n];
    const DG_FP y = v_y[n];
    v_x[n] = rx * x + ry * y;
    v_y[n] = sx * x + sy * y;
  }
}