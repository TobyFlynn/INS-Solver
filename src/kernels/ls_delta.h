inline void ls_delta(const DG_FP *alpha, const DG_FP *s, DG_FP *delta) {
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    // DG_FP del = 0.0;
    // if(fabs(s[i]) < *alpha) {
    //   del = 0.5 * ((1.0 / *alpha) + (1.0 / *alpha) * cos(PI * s[i] / *alpha));
    // }
    // delta[i] = del;
    // delta[i] = 0.5 * (PI / *alpha) * (1.0 / (cosh(PI * s[i] / *alpha) * cosh(PI * s[i] / *alpha)));

    DG_FP minDelta = 0.5 * (PI / *alpha) * (1.0 / (cosh(PI * 1.5) * cosh(PI * 1.5)));
    if(fabs(s[i]) < 1.5 * *alpha) {
      delta[i] = 0.5 * (PI / *alpha) * (1.0 / (cosh(PI * s[i] / *alpha) * cosh(PI * s[i] / *alpha))) - minDelta;
    } else {
      delta[i] = 0.0;
    }
  }
}
