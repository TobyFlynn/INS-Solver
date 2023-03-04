inline void init_surface_2d(const DG_FP *x, const DG_FP *y, DG_FP *s) {
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    // Level set: s is the distance from the interface
    s[i] = sqrt((x[i] * x[i]) / 2.0 + y[i] * y[i]) - 10.0;
  }
}
