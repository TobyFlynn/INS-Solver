inline void ins_2d_set_ls_vel(const DG_FP *t, const DG_FP *x, const DG_FP *y,
                              DG_FP *u, DG_FP *v) {
  const DG_FP _t = *t;
  for(int i = 0; i < DG_NP; i++) {
    ps2d_set_ls_vel(_t, x[i], y[i], u[i], v[i]);
  }
}
