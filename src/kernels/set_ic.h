inline void set_ic(const int *problem, const double *x, const double *y,
                   double *q0, double *q1) {
  const double PI = 3.141592653589793238463;
  const double R = 1.5;
  const double S = 13.5;
  for(int i = 0; i < DG_NP; i++) {
    double f = (1.0 - x[i] * x[i] - y[i] * y[i]) / (2.0 * R * R);
    q0[i] = (S * y[i] * exp(f)) / (2.0 * PI * R);
    q1[i] = (1.0 - ((S * x[i] * exp(f)) / (2.0 * PI * R)));
  }
}
