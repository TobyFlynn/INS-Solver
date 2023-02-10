inline void ins_set_art_vis(const double *f, const double *art_vis, double *vis, double *mm) {
  for(int i = 0; i < DG_NP; i++) {
    vis[i] = *art_vis + 1.0;
    mm[i] = *f;
  }
}