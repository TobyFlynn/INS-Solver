inline void sum_dg_np(DG_FP *res_g, const DG_FP *res) {
  for(int i = 0; i < DG_NP; i++) {
    *res_g += res[i];
  }
}
