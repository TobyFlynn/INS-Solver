inline void ins_3d_proj_1(DG_FP *tmp_alpha_0, DG_FP *tmp_alpha_1,
                          const DG_FP *r0, const DG_FP *r1, const DG_FP *r2,
                          const DG_FP *z0, const DG_FP *z1, const DG_FP *z2,
                          const DG_FP *p0, const DG_FP *p1, const DG_FP *p2,
                          const DG_FP *Ax0, const DG_FP *Ax1, const DG_FP *Ax2) {
  DG_FP a0_tmp = 0.0;
  DG_FP a1_tmp = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    a0_tmp += r0[i] * z0[i] + r1[i] * z1[i] + r2[i] * z2[i];
    a1_tmp += p0[i] * Ax0[i] + p1[i] * Ax1[i] + p2[i] * Ax2[i];
  }

  *tmp_alpha_0 += a0_tmp;
  *tmp_alpha_1 += a1_tmp;
}
