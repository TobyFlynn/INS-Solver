inline void ins_3d_bc_types(const DG_FP **nodes, int *type) {
  if(fp_equal(nodes[0][0], nodes[1][0]) && fp_equal(nodes[0][0], nodes[2][0]) && nodes[0][0] < -LW_INLET_LENGTH + 0.1) {
    *type = LW_INFLOW_BC;
  } else if(nodes[0][0] < 0.0 && nodes[1][0] < 0.0 && nodes[2][0] < 0.0) {
    *type = LW_NO_SLIP_WALL_BC;
  } else if(fp_equal(nodes[0][0], nodes[1][0]) && fp_equal(nodes[0][0], nodes[2][0]) && nodes[0][0] < LW_BLADE_START_X - 1e-1) {
    bool no_slip = false;
    for(int i = 0; i < 3; i++) {
      DG_FP dist = nodes[i][1] * nodes[i][1] + nodes[i][2] * nodes[i][2];
      if(dist < LW_INLET_NO_SLIP_RADIUS * LW_INLET_NO_SLIP_RADIUS)
        no_slip = true;
    }
    if(no_slip) {
      *type = LW_NO_SLIP_WALL_BC;
    } else {
      *type = LW_SLIP_WALL_BC;
    }
  } else if(fp_equal(nodes[0][0], nodes[1][0]) && fp_equal(nodes[0][0], nodes[2][0]) && nodes[0][0] > LW_LENGTH_SHORT - 0.1) {
    *type = LW_OUTFLOW_BC;
  } else if(nodes[0][1] * nodes[0][1] + nodes[0][2] * nodes[0][2] < LW_RADIUS * LW_RADIUS - 5e-1 &&
            nodes[1][1] * nodes[1][1] + nodes[1][2] * nodes[1][2] < LW_RADIUS * LW_RADIUS - 5e-1 &&
            nodes[2][1] * nodes[2][1] + nodes[2][2] * nodes[2][2] < LW_RADIUS * LW_RADIUS - 5e-1) {
    *type = LW_NO_SLIP_WALL_BC;
  } else {
    *type = LW_SLIP_WALL_BC;
  }
}
