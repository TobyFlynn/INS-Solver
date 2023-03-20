inline void ins_3d_bc_types(const DG_FP **nodes, int *type) {
  if(fp_equal(nodes[0][0], nodes[1][0]) && fp_equal(nodes[0][0], nodes[2][0]) && nodes[0][0] < -LW_INLET_LENGTH + 0.1) {
    *type = LW_INFLOW_BC;
  } else if(nodes[0][0] < 1e-8 && nodes[1][0] < 1e-8 && nodes[2][0] < 1e-8) {
    *type = LW_NO_SLIP_WALL_BC;
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
