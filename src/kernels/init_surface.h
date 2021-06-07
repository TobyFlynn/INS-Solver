inline void init_surface(const double *x, const double *y, double *s) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < 15; i++) {
    // Level set: s is the distance from the interface
    s[i] = sqrt((x[i] - 0.1) * (x[i] - 0.1) + (y[i] - 0.0) * (y[i] - 0.0)) - 0.15;
  }
}
