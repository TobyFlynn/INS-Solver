inline void diff_flux(const int *edgeNum, const bool *rev, const double **sJ,
                      const double **nx, const double **ny, const double **sigX,
                      const double **sigY, double **flux) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  // Do left element first
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

    double sigFX = (sigX[0][lInd] + sigX[1][rInd]) / 2.0;
    double sigFY = (sigY[0][lInd] + sigY[1][rInd]) / 2.0;
    // flux[0][exInd + i] += fscale[0][exInd + i] * (nx[0][exInd + i] * sigFX + ny[0][exInd + i] * sigFY - tau[0][edgeL] * (u[0][lInd] - u[1][rInd]));
    flux[0][lInd] += gaussW_g[i] * sJ[0][lInd] * (nx[0][lInd] * sigFX + ny[0][lInd] * sigFY);
  }

  for(int i = 0; i < 7; i++) {
    int lInd;
    int rInd = exIndR + i;
    if(reverse) {
      lInd = exIndL + 7 - i - 1;
    } else {
      lInd = exIndL + i;
    }

    double sigFX = (sigX[0][lInd] + sigX[1][rInd]) / 2.0;
    double sigFY = (sigY[0][lInd] + sigY[1][rInd]) / 2.0;
    // flux[1][exInd + i] += fscale[1][exInd + i] * (nx[1][exInd + i] * sigFX + ny[1][exInd + i] * sigFY - tau[1][edgeR] * (u[1][rInd] - u[0][lInd]));
    flux[1][rInd] += gaussW_g[i] * sJ[1][rInd] * (nx[1][rInd] * sigFX + ny[1][rInd] * sigFY);
  }
}
