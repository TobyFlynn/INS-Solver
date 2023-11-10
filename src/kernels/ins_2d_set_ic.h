inline void ins_2d_set_ic(const DG_FP *x, const DG_FP *y, DG_FP *u,
                          DG_FP *v, DG_FP *u1, DG_FP *v1) {
  for(int i = 0; i < DG_NP; i++) {
    ps2d_set_ic(x[i], y[i], u[i], v[i]);

    u1[i] = u[i];
    v1[i] = v[i];
  }
}
