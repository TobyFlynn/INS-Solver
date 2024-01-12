inline void mp_ins_3d_pr_2_oi(const DG_FP *fact, DG_FP *dx, DG_FP *dy, DG_FP *dz) {
  for(int i = 0; i < DG_CUB_3D_NP; i++) {
    dx[i] /= fact[i];
    dy[i] /= fact[i];
    dz[i] /= fact[i];
  }
}
