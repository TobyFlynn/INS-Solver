inline void diff_flux(const int *edgeNum, const bool *rev, const double **sJ,
                      const double **nx, const double **ny, const double **sigX,
                      const double **sigY, double **flux) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

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

    double sigFX = (sigX[0][lInd] + sigX[1][rInd]) / 2.0;
    double sigFY = (sigY[0][lInd] + sigY[1][rInd]) / 2.0;
    // flux[0][exInd + i] += fscale[0][exInd + i] * (nx[0][exInd + i] * sigFX + ny[0][exInd + i] * sigFY - tau[0][edgeL] * (u[0][lInd] - u[1][rInd]));
    flux[0][lInd] += gaussW_g[i] * sJ[0][lInd] * (nx[0][lInd] * sigFX + ny[0][lInd] * sigFY);
  }

  for(int i = 0; i < DG_GF_NP; i++) {
    int lInd;
    int rInd = exIndR + i;
    if(reverse) {
      lInd = exIndL + DG_GF_NP - i - 1;
    } else {
      lInd = exIndL + i;
    }

    double sigFX = (sigX[0][lInd] + sigX[1][rInd]) / 2.0;
    double sigFY = (sigY[0][lInd] + sigY[1][rInd]) / 2.0;
    // flux[1][exInd + i] += fscale[1][exInd + i] * (nx[1][exInd + i] * sigFX + ny[1][exInd + i] * sigFY - tau[1][edgeR] * (u[1][rInd] - u[0][lInd]));
    flux[1][rInd] += gaussW_g[i] * sJ[1][rInd] * (nx[1][rInd] * sigFX + ny[1][rInd] * sigFY);
  }
}
