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

  for(int m = 0; m < DG_GF_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind   = m * DG_NP + n;
      int indBC = *bedgeNum * DG_GF_NP + m;

      d[ind] = nx[indBC] * Dx[ind] + ny[indBC] * Dy[ind];
    }
  }
}
