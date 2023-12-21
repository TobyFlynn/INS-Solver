inline void mp_ins_2d_pr_2(const DG_FP *pr_factor, DG_FP *dpdx, DG_FP *dpdy) {
  for(int i = 0; i < DG_NP; i++) {
    dpdx[i] *= pr_factor[i];
    dpdy[i] *= pr_factor[i];
  }
}
