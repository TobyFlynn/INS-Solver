inline void ins_bc_types_2d(const DG_FP **nodes, int *type, int *pr_type,
                            int *vis_type) {
  // 0 is inflow
  // 1 is outflow
  // 2 is wall
  if(fabs(nodes[0][1] - nodes[1][1]) < 1e-8 && nodes[0][1] > 0.9) {
    *type = 1;
  } else {
    *type = 2;
  }

  if(*type == 1) {
    *pr_type = 0;
    *vis_type = 1;
  } else {
    *pr_type = 1;
    *vis_type = 0;
  }
}
