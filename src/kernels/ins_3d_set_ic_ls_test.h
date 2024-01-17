inline void ins_3d_set_ic_ls_test(const DG_FP *t, const DG_FP *x, const DG_FP *y, const DG_FP *z,
                                  DG_FP *u0, DG_FP *v0, DG_FP *w0, DG_FP *u1,
                                  DG_FP *v1, DG_FP *w1) {
  const DG_FP _t = *t;
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    const DG_FP time_component = cos(PI * _t / 8.0);
    u0[i] = sin(PI * x[i]) * sin(PI * x[i]) * sin(2.0 * PI * y[i]) * time_component;
    v0[i] = -sin(PI * y[i]) * sin(PI * y[i]) * sin(2.0 * PI * x[i]) * time_component;
    w0[i] = 0.0;
    u1[i] = u0[i];
    v1[i] = v0[i];
    w1[i] = w0[i];
  }
}
