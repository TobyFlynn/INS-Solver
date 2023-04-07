inline void ls_post_advec(const DG_FP *x, DG_FP *s) {
  for(int i = 0; i < DG_NP; i++) {
    if(x[i] < -LW_INLET_LENGTH + 0.5)
      s[i] = -1.0;
  }
}
