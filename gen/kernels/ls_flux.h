inline void ls_flux(const int *edgeNum, const bool *rev, const double **sJ,
                    const double **nx, const double **ny, const double **s,
                    double **dsldx, double **dsrdx, double **dsldy,
                    double **dsrdy) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  int exIndL = edgeL * 4;
  int exIndR = edgeR * 4;

  for(int i = 0; i < 4; i++) {
    int rInd;
    int lInd = exIndL + i;
    if(reverse) {
      rInd = exIndR + 4 - i - 1;
    } else {
      rInd = exIndR + i;
    }

    if(nx[0][lInd] >= 0.0) {
      dsldx[0][lInd] += gaussW_g[i] * sJ[0][lInd] * nx[0][lInd] * s[0][lInd];
      dsrdx[0][lInd] += gaussW_g[i] * sJ[0][lInd] * nx[0][lInd] * s[1][rInd];
    } else {
      dsldx[0][lInd] += gaussW_g[i] * sJ[0][lInd] * nx[0][lInd] * s[1][rInd];
      dsrdx[0][lInd] += gaussW_g[i] * sJ[0][lInd] * nx[0][lInd] * s[0][lInd];
    }

    if(ny[0][lInd] >= 0.0) {
      dsldy[0][lInd] += gaussW_g[i] * sJ[0][lInd] * ny[0][lInd] * s[0][lInd];
      dsrdy[0][lInd] += gaussW_g[i] * sJ[0][lInd] * ny[0][lInd] * s[1][rInd];
    } else {
      dsldy[0][lInd] += gaussW_g[i] * sJ[0][lInd] * ny[0][lInd] * s[1][rInd];
      dsrdy[0][lInd] += gaussW_g[i] * sJ[0][lInd] * ny[0][lInd] * s[0][lInd];
    }
  }

  for(int i = 0; i < 4; i++) {
    int lInd;
    int rInd = exIndR + i;
    if(reverse) {
      lInd = exIndL + 4 - i - 1;
    } else {
      lInd = exIndL + i;
    }

    if(nx[1][rInd] >= 0.0) {
      dsldx[1][rInd] += gaussW_g[i] * sJ[1][rInd] * nx[1][rInd] * s[1][rInd];
      dsrdx[1][rInd] += gaussW_g[i] * sJ[1][rInd] * nx[1][rInd] * s[0][lInd];
    } else {
      dsldx[1][rInd] += gaussW_g[i] * sJ[1][rInd] * nx[1][rInd] * s[0][lInd];
      dsrdx[1][rInd] += gaussW_g[i] * sJ[1][rInd] * nx[1][rInd] * s[1][rInd];
    }

    if(ny[1][rInd] >= 0.0) {
      dsldy[1][rInd] += gaussW_g[i] * sJ[1][rInd] * ny[1][rInd] * s[1][rInd];
      dsrdy[1][rInd] += gaussW_g[i] * sJ[1][rInd] * ny[1][rInd] * s[0][lInd];
    } else {
      dsldy[1][rInd] += gaussW_g[i] * sJ[1][rInd] * ny[1][rInd] * s[0][lInd];
      dsrdy[1][rInd] += gaussW_g[i] * sJ[1][rInd] * ny[1][rInd] * s[1][rInd];
    }
  }
}
