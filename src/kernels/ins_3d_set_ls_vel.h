inline void ins_3d_set_ls_vel(const DG_FP *t, const DG_FP *x, const DG_FP *y, const DG_FP *z,
                              DG_FP *u, DG_FP *v, DG_FP *w) {
  const DG_FP _t = *t;
  for(int i = 0; i < DG_NP; i++) {
    ps3d_set_ls_vel(_t, x[i], y[i], z[i], u[i], v[i], w[i]);
  }
}
