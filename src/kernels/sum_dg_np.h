inline void sum_dg_np(DG_FP *res_g, const DG_FP *res) {
  DG_FP tmp = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    tmp += res[i];
  }
  *res_g += tmp;
}
