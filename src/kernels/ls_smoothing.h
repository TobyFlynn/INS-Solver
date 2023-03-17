inline void ls_smoothing(const DG_FP *alpha, const DG_FP *s, DG_FP *smooth) {
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    // DG_FP step = 0.0;
    // if(s[i] > *alpha) {
    //   step = 1.0;
    // } else if(s[i] > -*alpha) {
    //   step = 0.5 * (1.0 + (s[i] / *alpha) + (1.0 / PI) * sin(PI * s[i] / *alpha));
    //   // step = 0.5 * (1.0 + (1.0 / PI) * sin(PI * s[i] / *alpha));
    // }
    // smooth[i] = step;
    smooth[i] = tanh(PI * s[i] / *alpha);

  }
}
