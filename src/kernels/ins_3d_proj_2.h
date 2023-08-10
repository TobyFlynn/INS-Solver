inline void ins_3d_proj_2(const DG_FP *_alpha, DG_FP *residual, DG_FP *tmp_beta_1,
                          const DG_FP *p0, const DG_FP *p1, const DG_FP *p2,
                          const DG_FP *Ax0, const DG_FP *Ax1, const DG_FP *Ax2,
                          DG_FP *z0, DG_FP *z1, DG_FP *z2, DG_FP *u, DG_FP *v,
                          DG_FP *w, DG_FP *r0, DG_FP *r1, DG_FP *r2) {
  const DG_FP alpha = *_alpha;
  DG_FP b_tmp = 0.0;
  DG_FP res_tmp = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    u[i] = u[i] + alpha * p0[i];
    v[i] = v[i] + alpha * p1[i];
    w[i] = w[i] + alpha * p2[i];

    b_tmp += r0[i] * z0[i] + r1[i] * z1[i] + r2[i] * z2[i];

    r0[i] = r0[i] - alpha * Ax0[i];
    r1[i] = r1[i] - alpha * Ax1[i];
    r2[i] = r2[i] - alpha * Ax2[i];

    res_tmp += r0[i] * r0[i] + r1[i] * r1[i] + r2[i] * r2[i];

    z0[i] = r0[i];
    z1[i] = r1[i];
    z2[i] = r2[i];
  }

  *tmp_beta_1 += b_tmp;
  *residual += res_tmp;
}
