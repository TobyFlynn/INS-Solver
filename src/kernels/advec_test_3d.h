inline void advec_test_3d(const DG_FP *x, const DG_FP *y, const DG_FP *z,
                          DG_FP *u, DG_FP *v, DG_FP *w, DG_FP *val) {
  for(int i = 0; i < DG_NP; i++) {
    u[i] = 1.0;
    v[i] = 1.0;
    w[i] = 1.0;
    val[i] = fmax(-(sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]) - (DG_FP)3.0), (DG_FP)0.0);
  }
}
