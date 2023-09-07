inline void ins_2d_proj_3(DG_FP *tmp_beta_0, const DG_FP *r0, const DG_FP *r1,
                          const DG_FP *z0, const DG_FP *z1) {
  DG_FP b0_tmp = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    b0_tmp += r0[i] * z0[i] + r1[i] * z1[i];
  }
  *tmp_beta_0 += b0_tmp;
}
