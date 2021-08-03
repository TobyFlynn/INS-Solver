inline void sigma_flux(const int *edgeNum, const bool *rev, const double **sJ,
                       const double **nx, const double **ny, const double **s,
                       double **sigFx, double **sigFy) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  // Do left element first
  int exIndL = edgeL * DG_GF_NP;
  int exIndR = edgeR * DG_GF_NP;

  for(int i = 0; i < DG_GF_NP; i++) {
    int rInd;
    int lInd = exIndL + i;
    if(reverse) {
      rInd = exIndR + DG_GF_NP - i - 1;
    } else {
      rInd = exIndR + i;
    }
    double flux = (s[0][lInd] + s[1][rInd]) / 2.0;
    sigFx[0][lInd] += gaussW_g[i] * sJ[0][lInd] * nx[0][lInd] * flux;
    sigFy[0][lInd] += gaussW_g[i] * sJ[0][lInd] * ny[0][lInd] * flux;
  }

  for(int i = 0; i < DG_GF_NP; i++) {
    int lInd;
    int rInd = exIndR + i;
    if(reverse) {
      lInd = exIndL + DG_GF_NP - i - 1;
    } else {
      lInd = exIndL + i;
    }
    double flux = (s[0][lInd] + s[1][rInd]) / 2.0;
    sigFx[1][rInd] += gaussW_g[i] * sJ[1][rInd] * nx[1][rInd] * flux;
    sigFy[1][rInd] += gaussW_g[i] * sJ[1][rInd] * ny[1][rInd] * flux;
  }
}
