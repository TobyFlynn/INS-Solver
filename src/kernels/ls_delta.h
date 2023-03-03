inline void ls_delta(const DG_FP *alpha, const DG_FP *s, DG_FP *delta) {
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    if(fabs(s[i]) < *alpha)
      delta[i] = 0.5 * (*alpha) * (1.0 + cos(PI * s[i] / *alpha));
    else
      delta[i] = 0.0;
  }
}
