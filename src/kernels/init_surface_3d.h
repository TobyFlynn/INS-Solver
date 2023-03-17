inline void init_surface_3d(const DG_FP *x, const DG_FP *y, const DG_FP *z,
                            DG_FP *s) {
  const DG_FP PI = 3.141592653589793238463;

  for(int i = 0; i < DG_NP; i++) {
    // s[i] = sqrt((x[i] - 3.0) * (x[i] - 3.0) + (y[i] - 3.0) * (y[i] - 3.0) + (z[i] - 3.0) * (z[i] - 3.0)) - 1.5;
    // s[i] = x[i] - 1.0;
    // s[i] = sqrt((x[i] - 1.0) * (x[i] - 1.0) + (y[i] - 0.5) * (y[i] - 0.5) + (z[i] - 0.5) * (z[i] - 0.5)) - 0.25;
    // s[i] = sqrt((x[i] - 0.1) * (x[i] - 0.1) + y[i] * y[i] + z[i] * z[i]) - 0.05;
    s[i] = sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]) - LW_INLET_RADIUS;
  }
}
