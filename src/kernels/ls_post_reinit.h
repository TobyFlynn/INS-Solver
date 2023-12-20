inline void ls_post_reinit(DG_FP *s) {
  for(int i = 0; i < DG_NP; i++) {
    s[i] = fmax(fmin(LS_CAP, s[i]), -LS_CAP);
  }
}
