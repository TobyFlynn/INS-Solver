inline void init_gauss_grad4(const int *bedgeNum,
                             const double *nx, const double *ny,
                             const double *Dx0, const double *Dy0,
                             const double *Dx1, const double *Dy1,
                             const double *Dx2, const double *Dy2, double *d) {
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

  for(int j = 0; j < 10; j++) {
    for(int i = 0; i < 6; i++) {
      int ind   = j * 6 + i;
      int indBC = *bedgeNum * 6 + i;

      d[ind] = nx[indBC] * Dx[ind] + ny[indBC] * Dy[ind];
    }
  }
}
