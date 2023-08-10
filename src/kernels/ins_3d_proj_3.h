inline void ins_3d_proj_3(DG_FP *tmp_beta_0, const DG_FP *r0, const DG_FP *r1,
                          const DG_FP *r2, const DG_FP *z0, const DG_FP *z1,
                          const DG_FP *z2) {
  DG_FP b0_tmp = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    b0_tmp += r0[i] * z0[i] + r1[i] * z1[i] + r2[i] * z2[i];
  }
  *tmp_beta_0 += b0_tmp;
}
