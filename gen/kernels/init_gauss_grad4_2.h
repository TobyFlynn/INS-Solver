inline void init_gauss_grad4_2(const int *bedgeNum,
                             const double *nx, const double *ny,
                             const double *Dx0, const double *Dy0,
                             const double *Dx1, const double *Dy1,
                             const double *Dx2, const double *Dy2,
                             const double *fact, double *d) {
  const double *Dx, *Dy;

  if(*bedgeNum == 0) {
    Dx = Dx0;
    Dy = Dy0;
  } else if(*bedgeNum == 1) {
    Dx = Dx1;
    Dy = Dy1;
  } else {
    Dx = Dx2;
    Dy = Dy2;
  }

  for(int m = 0; m < 6; m++) {
    for(int n = 0; n < 10; n++) {
      int ind   = m * 10 + n;
      int indBC = *bedgeNum * 6 + m;

      d[ind] = nx[indBC] * fact[indBC] * Dx[ind] + ny[indBC] * fact[indBC] * Dy[ind];
    }
  }
}
