inline void init_surface(const double *x, const double *y, double *s) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    // Level set: s is the distance from the interface
    s[i] = sqrt(x[i] * x[i] + y[i] * y[i]) - 3.0;
  }
}
