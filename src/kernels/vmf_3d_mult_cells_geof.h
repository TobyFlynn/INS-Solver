inline void vmf_3d_mult_cells_geof(const DG_FP *geof, DG_FP *u_x, DG_FP *u_y, 
                          DG_FP *u_z, DG_FP *v_x, DG_FP *v_y, DG_FP *v_z, 
                          DG_FP *w_x, DG_FP *w_y, DG_FP *w_z) {
  for(int n = 0; n < DG_NP; n++) {
    const DG_FP ux = u_x[n];
    const DG_FP uy = u_y[n];
    const DG_FP uz = u_z[n];
    u_x[n] = geof[RX_IND] * ux + geof[RY_IND] * uy + geof[RZ_IND] * uz;
    u_y[n] = geof[SX_IND] * ux + geof[SY_IND] * uy + geof[SZ_IND] * uz;
    u_z[n] = geof[TX_IND] * ux + geof[TY_IND] * uy + geof[TZ_IND] * uz;
    const DG_FP vx = v_x[n];
    const DG_FP vy = v_y[n];
    const DG_FP vz = v_z[n];
    v_x[n] = geof[RX_IND] * vx + geof[RY_IND] * vy + geof[RZ_IND] * vz;
    v_y[n] = geof[SX_IND] * vx + geof[SY_IND] * vy + geof[SZ_IND] * vz;
    v_z[n] = geof[TX_IND] * vx + geof[TY_IND] * vy + geof[TZ_IND] * vz;
    const DG_FP wx = w_x[n];
    const DG_FP wy = w_y[n];
    const DG_FP wz = w_z[n];
    w_x[n] = geof[RX_IND] * wx + geof[RY_IND] * wy + geof[RZ_IND] * wz;
    w_y[n] = geof[SX_IND] * wx + geof[SY_IND] * wy + geof[SZ_IND] * wz;
    w_z[n] = geof[TX_IND] * wx + geof[TY_IND] * wy + geof[TZ_IND] * wz;
  }
}
