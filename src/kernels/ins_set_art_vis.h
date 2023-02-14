inline void ins_set_art_vis(const DG_FP *f, const DG_FP *art_vis, DG_FP *vis, DG_FP *mm) {
  for(int i = 0; i < DG_NP; i++) {
    vis[i] = *art_vis + 1.0;
    mm[i] = *f;
  }
}