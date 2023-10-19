inline void ins_2d_pr_bc(const int *bc_type, const int *type, DG_FP *bcs) {
  for(int i = 0; i < DG_NPF; i++) {
    bcs[i] = 0.0;
  }
}
