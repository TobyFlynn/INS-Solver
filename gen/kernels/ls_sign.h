inline void ls_sign(const double *a, const double *s, const double *dsdx,
                    const double *dsdy, double *sign) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < 6; i++) {
    double alpha = *a * sqrt(dsdx[i] * dsdx[i] + dsdy[i] * dsdy[i]);
    sign[i] = tanh(PI * s[i] / alpha);
  }
}
