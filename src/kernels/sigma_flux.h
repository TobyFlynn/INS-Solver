inline void sigma_flux(const int **p, const int *edgeNum, const bool *rev,
                       const double **sJ, const double **nx, const double **ny,
                       const double **s, const double **vis, double **sigFx,
                       double **sigFy) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  // Get constants
  const double *gaussW = &gaussW_g[(p[0][0] - 1) * DG_GF_NP];

  // Do left element first
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
  }

  for(int i = 0; i < DG_GF_NP; i++) {
    int rInd;
    int lInd = exIndL + i;
    if(reverse) {
      rInd = exIndR + DG_GF_NP - i - 1;
    } else {
      rInd = exIndR + i;
    }
    double flux = wL * s[0][lInd] + wR * s[1][rInd];
    flux *= kL;
    sigFx[0][lInd] += gaussW[i] * sJ[0][lInd] * nx[0][lInd] * flux;
    sigFy[0][lInd] += gaussW[i] * sJ[0][lInd] * ny[0][lInd] * flux;
  }

  for(int i = 0; i < DG_GF_NP; i++) {
    int lInd;
    int rInd = exIndR + i;
    if(reverse) {
      lInd = exIndL + DG_GF_NP - i - 1;
    } else {
      lInd = exIndL + i;
    }
    double flux = wL * s[0][lInd] + wR * s[1][rInd];
    flux *= kR;
    sigFx[1][rInd] += gaussW[i] * sJ[1][rInd] * nx[1][rInd] * flux;
    sigFy[1][rInd] += gaussW[i] * sJ[1][rInd] * ny[1][rInd] * flux;
  }
}
