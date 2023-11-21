inline void ins_2d_set_vis_y_bc_type(const int *type, int *vis_type) {
  if(*type == BC_TYPE_NO_SLIP) {
    *vis_type = BC_DIRICHLET;
  } else if(*type == BC_TYPE_SLIP_X) {
    *vis_type = BC_NEUMANN;
  } else if(*type == BC_TYPE_SLIP_Y) {
    *vis_type = BC_DIRICHLET;
  } else if(*type == BC_TYPE_NATURAL_OUTFLOW) {
    *vis_type = BC_NEUMANN;
  } else {
    ps2d_custom_bc_get_vis_type(*type, *vis_type);
  }
}
