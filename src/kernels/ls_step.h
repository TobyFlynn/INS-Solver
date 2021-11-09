inline void ls_step(const double *alpha, const double *s, double *step) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    step[i] = tanh(PI * s[i] / *alpha);
  }
}
