inline void mp_ins_2d_pr_oi_1(DG_FP *out) {
  for(int i = 0; i < DG_NUM_FACES * DG_CUB_SURF_2D_NP; i++) {
    out[i] = 1.0 / out[i];
  }
}
