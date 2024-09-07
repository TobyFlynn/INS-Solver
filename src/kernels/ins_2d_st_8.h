inline void ins_2d_st_8(const DG_FP *h, const DG_FP *s, DG_FP *curv) {
  for(int i = 0; i < DG_NP; i++) {
    if(fabs(curv[i]) > 1e-6 && fabs((1.0 / curv[i]) - s[i]) > 1e-6)
      curv[i] = 1.0 / ((1.0 / curv[i]) - s[i]);
    if(fabs(curv[i]) > 1.0 / *h)
      curv[i] = curv[i] > 0.0 ? 1.0 / *h : -1.0 / *h;
  }
}
