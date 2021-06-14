inline void poisson_mf_edges(const int *edgeNum, const bool *rev,
                             const double **sJ, const double **nx,
                             const double **ny, const double **tau,
                             const double **u, const double **dudx,
                             const double **dudy, double **fluxX,
                             double **fluxY, double **flux) {
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  int exIndL = 0;
  if(edgeL == 1) exIndL = 7;
  else if(edgeL == 2) exIndL = 2 * 7;

  int exIndR = 0;
  if(edgeR == 1) exIndR = 7;
  else if(edgeR == 2) exIndR = 2 * 7;

  for(int i = 0; i < 7; i++) {
    int rInd;
    int lInd = exIndL + i;
    if(reverse) {
      rInd = exIndR + 7 - i - 1;
    } else {
      rInd = exIndR + i;
    }

    double tmp = (u[0][lInd] + u[1][rInd]) / 2.0;
    tmp *= gaussW_g[i] * sJ[0][lInd];
    fluxX[0][lInd] += nx[0][lInd] * tmp;
    fluxY[0][lInd] += ny[0][lInd] * tmp;
    tmp = nx[0][lInd] * ((dudx[0][lInd] + dudx[1][rInd]) / 2.0);
    tmp += ny[0][lInd] * ((dudy[0][lInd] + dudy[1][rInd]) / 2.0);
    tmp -= tau[0][edgeL] * (u[0][lInd] - u[1][rInd]) / 2.0;
    tmp *= gaussW_g[i] * sJ[0][lInd];
    flux[0][lInd] += tmp;
  }

  for(int i = 0; i < 7; i++) {
    int lInd;
    int rInd = exIndR + i;
    if(reverse) {
      lInd = exIndL + 7 - i - 1;
    } else {
      lInd = exIndL + i;
    }

    double tmp = (u[0][lInd] + u[1][rInd]) / 2.0;
    tmp *= gaussW_g[i] * sJ[1][rInd];
    fluxX[1][rInd] += nx[1][rInd] * tmp;
    fluxY[1][rInd] += ny[1][rInd] * tmp;
    tmp = nx[1][rInd] * ((dudx[0][lInd] + dudx[1][rInd]) / 2.0);
    tmp += ny[1][rInd] * ((dudy[0][lInd] + dudy[1][rInd]) / 2.0);
    tmp -= tau[1][edgeR] * (u[1][rInd] - u[0][lInd]) / 2.0;
    tmp *= gaussW_g[i] * sJ[1][rInd];
    flux[1][rInd] += tmp;
  }
}
