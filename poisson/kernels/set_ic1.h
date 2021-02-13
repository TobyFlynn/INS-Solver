inline void set_ic1(double *uD, double *qN) {
  for(int i = 0; i < 15; i++) {
    uD[i] = 0.0;
    qN[i] = 0.0;
  }
}
