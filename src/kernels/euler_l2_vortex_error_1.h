inline void euler_l2_vortex_error_1(double *res_g, const double *res) {
  for(int i = 0; i < DG_NP; i++) {
    *res_g += res[i];
  }
}
