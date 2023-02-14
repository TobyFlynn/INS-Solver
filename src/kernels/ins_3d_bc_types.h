inline void ins_3d_bc_types(const DG_FP **nodes, int *type) {
  // 0 is inflow
  // 1 is outflow
  // 2 is wall
  const DG_FP inlet_radius_2 = 0.05 * 0.05 + 1e-8;
  const DG_FP inlet_radius_2_ = 0.2 * 0.2 + 1e-8;

  if(fabs(nodes[0][0] - nodes[1][0]) < 1e-8 && fabs(nodes[0][0] - nodes[2][0]) < 1e-8 && nodes[0][0] < 0.1) {
    DG_FP tmp0 = nodes[0][0] * nodes[0][0] + nodes[0][1] * nodes[0][1] + nodes[0][2] * nodes[0][2];
    DG_FP tmp1 = nodes[1][0] * nodes[1][0] + nodes[1][1] * nodes[1][1] + nodes[1][2] * nodes[1][2];
    DG_FP tmp2 = nodes[2][0] * nodes[2][0] + nodes[2][1] * nodes[2][1] + nodes[2][2] * nodes[2][2];
    if(tmp0 < inlet_radius_2 && tmp1 < inlet_radius_2 && tmp2 < inlet_radius_2) {
      *type = 0;
    } else if(tmp0 < inlet_radius_2_ && tmp1 < inlet_radius_2_ && tmp2 < inlet_radius_2_) {
      *type = 5;
    } else {
      *type = 2;
    }
  } else if(fabs(nodes[0][0] - nodes[1][0]) < 1e-8 && fabs(nodes[0][0] - nodes[2][0]) < 1e-8 && nodes[0][0] > 0.9) {
    *type = 1;
  } else {
    *type = 2;
  }
}
