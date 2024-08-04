inline void cg_update_ans(const DG_FP *alpha, const DG_FP *p0, const DG_FP *p1,
                          DG_FP *u, DG_FP *v) {
  for(int i = 0; i < DG_NP; i++) {
    u[i] += *alpha * p0[i];
    v[i] += *alpha * p1[i];
  }
}