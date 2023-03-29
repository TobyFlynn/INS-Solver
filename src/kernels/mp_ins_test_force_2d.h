inline void mp_ins_test_force_2d(const double *x, const double *y,
                                 double *sf_x, double *sf_y) {
  const DG_FP a = 1.0;
  const DG_FP b = 0.0;
  const DG_FP c = 0.1;
  for(int i = 0; i < DG_NP; i++) {
    sf_x[i] = 0.0;
    sf_y[i] = a * exp((-(x[i] - b) * (x[i] - b)) / (2.0 * c * c));
  }
}
