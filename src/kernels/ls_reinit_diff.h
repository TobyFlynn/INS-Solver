inline void ls_reinit_diff(const DG_FP *x, const DG_FP *y, DG_FP *s) {
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    // Level set: s is the distance from the interface
    s[i] = fabs(s[i] - (sqrt((x[i] - 0.5) * (x[i] - 0.5) + (y[i] - 0.75) * (y[i] - 0.75)) - 0.1));
  }
}
