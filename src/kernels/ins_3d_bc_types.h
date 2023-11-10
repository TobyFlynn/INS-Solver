inline void ins_3d_bc_types(const DG_FP **nodes, int *type) {
  ps3d_set_boundary_type(nodes[0][0], nodes[0][1], nodes[0][2], nodes[1][0],
        nodes[1][1], nodes[1][2], nodes[2][0], nodes[2][1], nodes[2][2], *type);
}
