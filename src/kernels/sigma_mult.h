inline void sigma_mult(const double *eps, double *sigx, double *sigy,
                       double *diffF) {
  for(int i = 0; i < DG_NP; i++) {
    sigx[i] *= *eps;
    sigy[i] *= *eps;
  }

  for(int i = 0; i < DG_G_NP; i++) {
    diffF[i] = 0.0;
  }
}
