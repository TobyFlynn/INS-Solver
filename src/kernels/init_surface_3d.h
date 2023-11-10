inline void init_surface_3d(const DG_FP *x, const DG_FP *y, const DG_FP *z,
                            DG_FP *s) {
  for(int i = 0; i < DG_NP; i++) {
    ps3d_set_surface(x[i], y[i], z[i], s[i]);
  }
}
