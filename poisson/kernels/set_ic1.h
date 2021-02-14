inline void set_ic1(double *uD, double *qN, double *rhs, double *tau,
                    double *ex1, double *ex2) {
  for(int i = 0; i < 15; i++) {
    uD[i] = 0.0;
    qN[i] = 0.0;
    rhs[i] = 0.0;
    tau[i] = 0.0;
    ex1[i] = 0.0;
    ex2[i] = 0.0;
  }
}
