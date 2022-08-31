inline void set_ic(const int *problem, const double *x, const double *y,
                   double *q0, double *q1) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    q0[i] = 1.0;
    q1[i] = 0.0;
  }
}
