inline void mp_ins_pr_fsv(const DG_FP *factor, DG_FP *dpdy) {
  const DG_FP fact = *factor;
  for(int i = 0; i < DG_NP; i++) {
    dpdy[i] += fact;
  }
}
