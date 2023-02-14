inline void euler_2d_l2_vortex_error_1(DG_FP *res_g, const DG_FP *res) {
  for(int i = 0; i < DG_NP; i++) {
    *res_g += res[i];
  }
}