inline void set_ic1(double *uD, double *qN, double *rhs) {
  for(int i = 0; i < 15; i++) {
    uD[i] = 0.0;
    qN[i] = 0.0;
    rhs[i] = 0.0;
  }
}
