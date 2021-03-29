inline void diff(const double *x, double *res) {
  for(int i = 0; i < 15; i++) {
    res[i] = x[i] - res[i];
    // res[i] = x[i];
  }
}
