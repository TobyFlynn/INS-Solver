inline void mp_ins_3d_pr_fact_oi_1(DG_FP *out) {
  for(int i = 0; i < DG_NUM_FACES * DG_CUB_SURF_3D_NP; i++) {
    out[i] = 1.0 / out[i];
  }
}
