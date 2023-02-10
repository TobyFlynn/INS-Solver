inline void ins_3d_pr_2(const int *bc_type, int *type, double *bcs) {
  if(*bc_type == 1)
    *type = 0;
  else
    *type = 1;
  for(int i = 0; i < DG_NPF; i++) {
    bcs[i] = 0.0;
  }
}
