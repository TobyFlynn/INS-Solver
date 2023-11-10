inline void ins_bc_types(const DG_FP **nodes, int *type, int *pr_type,
                         int *vis_type) {
  ps2d_set_boundary_type(nodes[0][0], nodes[0][1], nodes[1][0], nodes[1][1], *type);

  if(*type == BC_TYPE_NO_SLIP) {
    *pr_type = BC_NEUMANN;
    *vis_type = BC_DIRICHLET;
  } else if(*type == BC_TYPE_SLIP) {
    *pr_type = BC_NEUMANN;
    *vis_type = BC_DIRICHLET;
  } else if(*type == BC_TYPE_NATURAL_OUTFLOW) {
    *pr_type = BC_DIRICHLET;
    *vis_type = BC_NEUMANN;
  } else {
    ps2d_custom_bc_get_pr_type(*type, *pr_type);
    ps2d_custom_bc_get_vis_type(*type, *vis_type);
  }
}
