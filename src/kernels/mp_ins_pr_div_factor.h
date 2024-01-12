inline void mp_ins_pr_div_factor(const DG_FP *_div_factor, DG_FP *div) {
  const DG_FP div_factor = *_div_factor;
  for(int i = 0; i < DG_NP; i++) {
    div[i] *= div_factor;
  }
}
