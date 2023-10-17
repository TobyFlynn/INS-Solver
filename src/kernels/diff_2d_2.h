inline void diff_2d_2(const DG_FP *vis, DG_FP *val_x, DG_FP *val_y) {
  for(int i = 0; i < DG_NP; i++) {
    val_x[i] *= vis[i];
    val_y[i] *= vis[i];
  }
}