inline void sigma_mult(const double *eps, double *sigx, double *sigy,
                       double *fx, double *fy, double *diffF) {
  for(int i = 0; i < 10; i++) {
    sigx[i] *= *eps;
    sigy[i] *= *eps;
  }

  for(int i = 0; i < 18; i++) {
    fx[i] = 0.0;
    fy[i] = 0.0;
    diffF[i] = 0.0;
  }
}
