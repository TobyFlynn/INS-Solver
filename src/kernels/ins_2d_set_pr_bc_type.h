inline void ins_2d_set_pr_bc_type(const int *type, int *pr_type) {
  if(*type == BC_TYPE_NO_SLIP) {
    *pr_type = BC_NEUMANN;
  } else if(*type == BC_TYPE_SLIP) {
    *pr_type = BC_NEUMANN;
  } else if(*type == BC_TYPE_NATURAL_OUTFLOW) {
    *pr_type = BC_DIRICHLET;
  } else {
    ps2d_custom_bc_get_pr_type(*type, *pr_type);
  }
}
