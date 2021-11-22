inline void diff_flux(const int *edgeNum, const bool *rev, const double **sJ,
                      const double **nx, const double **ny, const double **s,
                      const double **vis, const double **sigX,
                      const double **sigY, double **flux) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  int exIndL = edgeL * DG_GF_NP;
  int exIndR = edgeR * DG_GF_NP;

  double kL = sqrt(vis[0][0]);
  double kR = sqrt(vis[1][0]);
  double lamdaL = kL;
  double lamdaR = kR;
  double wL, wR;

  if(lamdaL < 1e-12 && lamdaR < 1e-12) {
    wL = 0.5;
    wR = 0.5;
  } else {
    double lamdaAvg = (lamdaL + lamdaR) / 2.0;
    wL = lamdaL / (2.0 * lamdaAvg);
    wR = lamdaR / (2.0 * lamdaAvg);

    wL = 1.0 - wL;
    wR = 1.0 - wR;
  }

  for(int i = 0; i < DG_GF_NP; i++) {
    int rInd;
    int lInd = exIndL + i;
    if(reverse) {
      rInd = exIndR + DG_GF_NP - i - 1;
    } else {
      rInd = exIndR + i;
    }

    double sigFX = wL * sigX[0][lInd] + wR * sigX[1][rInd];
    double sigFY = wL * sigY[0][lInd] + wR * sigY[1][rInd];
    // This pen term needs checking
    double pen   = fmin(lamdaL, lamdaR) * fmin(lamdaL, lamdaR) * (s[0][lInd] - s[1][rInd]);
    // flux[0][exInd + i] += fscale[0][exInd + i] * (nx[0][exInd + i] * sigFX + ny[0][exInd + i] * sigFY - tau[0][edgeL] * (u[0][lInd] - u[1][rInd]));
    flux[0][lInd] += gaussW_g[i] * sJ[0][lInd] * (nx[0][lInd] * sigFX + ny[0][lInd] * sigFY - pen);
  }

  for(int i = 0; i < DG_GF_NP; i++) {
    int lInd;
    int rInd = exIndR + i;
    if(reverse) {
      lInd = exIndL + DG_GF_NP - i - 1;
    } else {
      lInd = exIndL + i;
    }

    double sigFX = wL * sigX[0][lInd] + wR * sigX[1][rInd];
    double sigFY = wL * sigY[0][lInd] + wR * sigY[1][rInd];
    // This pen term needs checking
    double pen   = fmin(lamdaL, lamdaR) * fmin(lamdaL, lamdaR) * (s[1][rInd] - s[0][lInd]);
    // flux[1][exInd + i] += fscale[1][exInd + i] * (nx[1][exInd + i] * sigFX + ny[1][exInd + i] * sigFY - tau[1][edgeR] * (u[1][rInd] - u[0][lInd]));
    flux[1][rInd] += gaussW_g[i] * sJ[1][rInd] * (nx[1][rInd] * sigFX + ny[1][rInd] * sigFY - pen);
  }
}
