inline void ins_3d_bc_types(const DG_FP **nodes, int *type) {
  if(fp_equal(nodes[0][0], nodes[1][0]) && fp_equal(nodes[0][0], nodes[2][0]) && nodes[0][0] < 0.1) {
    bool inlet = true;
    for(int i = 0; i < 3; i++) {
      DG_FP dist = nodes[i][1] * nodes[i][1] + nodes[i][2] * nodes[i][2];
      if(dist > LW_INLET_RADIUS * LW_INLET_RADIUS + 1e-8)
        inlet = false;
    }
    if(inlet) {
      *type = LW_INFLOW_BC;
    } else {
      *type = LW_SLIP_WALL_BC;
    }
  } else if(fp_equal(nodes[0][0], nodes[1][0]) && fp_equal(nodes[0][0], nodes[2][0]) && nodes[0][0] > LW_LENGTH - 0.1) {
    *type = LW_OUTFLOW_BC;
  } else {
    *type = LW_SLIP_WALL_BC;
  }
}
