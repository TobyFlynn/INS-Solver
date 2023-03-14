inline void ins_2d_set_ic(const DG_FP *x, const DG_FP *y, DG_FP *u,
                          DG_FP *v, DG_FP *u1, DG_FP *v1) {
  for(int i = 0; i < DG_NP; i++) {
    u[i] = 0.0;
    v[i] = 0.0;
    u1[i] = u[i];
    v1[i] = v[i];
  }
}
