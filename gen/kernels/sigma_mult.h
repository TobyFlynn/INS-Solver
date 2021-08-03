inline void sigma_mult(const double *eps, double *sigx, double *sigy,
                       double *diffF) {
  for(int i = 0; i < 10; i++) {
    sigx[i] *= *eps;
    sigy[i] *= *eps;
  }

  for(int i = 0; i < 18; i++) {
    diffF[i] = 0.0;
  }
}
