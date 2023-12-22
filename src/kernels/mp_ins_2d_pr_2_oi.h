inline void mp_ins_2d_pr_2_oi(const DG_FP *rho, DG_FP *dx, DG_FP *dy) {
  for(int i = 0; i < DG_CUB_2D_NP; i++) {
    dx[i] /= rho[i];
    dy[i] /= rho[i];
  }
}
