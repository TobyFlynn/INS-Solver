inline void cg3_update_ans(const DG_FP *alpha, const DG_FP *p0, const DG_FP *p1,
                           const DG_FP *p2, DG_FP *u, DG_FP *v, DG_FP *w) {
  for(int i = 0; i < DG_NP; i++) {
    u[i] += *alpha * p0[i];
    v[i] += *alpha * p1[i];
    w[i] += *alpha * p2[i];
  }
}