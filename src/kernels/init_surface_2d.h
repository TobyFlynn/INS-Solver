inline void init_surface_2d(const DG_FP *x, const DG_FP *y, DG_FP *s) {
  for(int i = 0; i < DG_NP; i++) {
    // Level set: s is the distance from the interface
    ps2d_set_surface(x[i], y[i], s[i]);
  }
}
