inline void ins_2d_st_8(const DG_FP *alpha, const DG_FP *s, DG_FP *curv) {
  const DG_FP _alpha = *alpha;
  for(int i = 0; i < DG_NP; i++) {
    if(fabs(s[i]) < _alpha)
      curv[i] = 1.0 / ((1.0 / curv[i]) - s[i]);
  }
}