inline void ins_2d_st_art_vis(const DG_FP *alpha_, const DG_FP *s, DG_FP *out) {
  const DG_FP alpha = *alpha_ * 2.0;
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    out[i] = fmax(1e-6, fmin(1.0 , 2.0 * (alpha - fabs(s[i])) / alpha));
  }
}