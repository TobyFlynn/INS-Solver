inline void cg3_compute_beta_pre_flex(DG_FP *beta, const DG_FP *r0, const DG_FP *r1, const DG_FP *r2,
                                      const DG_FP *z0, const DG_FP *z1, const DG_FP *z2,
                                      const DG_FP *z0_1, const DG_FP *z1_1, const DG_FP *z2_1) {
  DG_FP tmp = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    tmp += r0[i] * (z0[i] - z0_1[i]) + r1[i] * (z1[i] - z1_1[i]) + r2[i] * (z2[i] - z2_1[i]);
  }
  *beta += tmp;
}