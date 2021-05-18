inline void set_ic(const int *problem, const double *x, const double *y,
                   double *q0, double *q1) {
  const double PI = 3.141592653589793238463;
  if(*problem == 0) {
    for(int i = 0; i < 15; i++) {
      q0[i] = 0.0;
      q1[i] = 0.0;
    }
  } else {
    for(int i = 0; i < 15; i++) {
      q0[i] = -sin(2.0 * PI * y[i]) * exp(-nu * 4.0 * PI * PI * 0.0);
      q1[i] = sin(2.0 * PI * x[i]) * exp(-nu * 4.0 * PI * PI * 0.0);
    }
  }
}
