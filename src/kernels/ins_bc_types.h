inline void ins_bc_types(const DG_FP **nodes, int *type) {
  ps2d_set_boundary_type(nodes[0][0], nodes[0][1], nodes[1][0], nodes[1][1], *type);
}
