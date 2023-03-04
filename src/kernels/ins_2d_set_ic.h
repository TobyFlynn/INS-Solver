inline void ins_2d_set_ic(const DG_FP *x, const DG_FP *y, DG_FP *u,
                          DG_FP *v, DG_FP *u1, DG_FP *v1) {
  const DG_FP PI = 3.141592653589793238463;
  const DG_FP R = 1.5;
  const DG_FP S = 13.5;
  for(int i = 0; i < DG_NP; i++) {
    DG_FP f = (1.0 - x[i] * x[i] - y[i] * y[i]) / (2.0 * R * R);
    // u[i] = (S * y[i] * exp(f)) / (2.0 * PI * R);
    // v[i] = (1.0 - ((S * x[i] * exp(f)) / (2.0 * PI * R)));
    u[i] = 1.0;
    v[i] = 0.0;
    u1[i] = u[i];
    v1[i] = v[i];
  }
}
