inline void set_ic(const int *problem, const double *x, const double *y,
                   double *q0, double *q1) {
  const double PI = 3.141592653589793238463;
  if(*problem == 0) {
    for(int i = 0; i < DG_NP; i++) {
      // q0[i] = 0.1;
      q0[i] = sin((PI * 1.0) / 8.0) * 4.0 * y[i] * (1.0 - y[i]);
      q1[i] = 0.0;
    }
  } else if(*problem == 1) {
    for(int i = 0; i < DG_NP; i++) {
      // q0[i] = 0.0;
      // q1[i] = 0.0;
      q0[i] = sin(PI * x[i]) * sin(PI * x[i]) * sin(2.0 * PI * y[i]);
      q1[i] = -sin(2.0 * PI * x[i]) * sin(PI * y[i]) * sin(PI * y[i]);
    }
  } else {
    for(int i = 0; i < DG_NP; i++) {
      q0[i] = -sin(2.0 * PI * y[i]) * exp(-nu * 4.0 * PI * PI * 0.0);
      q1[i] = sin(2.0 * PI * x[i]) * exp(-nu * 4.0 * PI * PI * 0.0);
    }
  }
}
