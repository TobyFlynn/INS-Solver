inline void setup_poisson(double *tau, double *ex1, double *ex2) {
  for(int i = 0; i < 15; i++) {
    tau[i] = 0.0;
    ex1[i] = 0.0;
    ex2[i] = 0.0;
  }
}
