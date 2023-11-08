inline void ins_3d_set_ic(const DG_FP *x, const DG_FP *y, const DG_FP *z,
                          DG_FP *u0, DG_FP *v0, DG_FP *w0, DG_FP *u1,
                          DG_FP *v1, DG_FP *w1) {
  for(int i = 0; i < DG_NP; i++) {
    ps3d_set_ic(x[i], y[i], z[i], u0[i], v0[i], w0[i]);

    u1[i] = u0[i];
    v1[i] = v0[i];
    w1[i] = w0[i];
  }
}
