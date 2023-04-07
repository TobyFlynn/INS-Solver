inline void ls_post_reinit(DG_FP *s) {
  for(int i = 0; i < DG_NP; i++) {
    s[i] = fmax(fmin(1.0, s[i]), -1.0);
  }
}
