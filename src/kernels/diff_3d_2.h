inline void diff_3d_2(const DG_FP *vis, DG_FP *val_x, DG_FP *val_y, DG_FP *val_z) {
  for(int i = 0; i < DG_NP; i++) {
    const DG_FP _vis = vis[i];
    val_x[i] *= _vis;
    val_y[i] *= _vis;
    val_z[i] *= _vis;
  }
}