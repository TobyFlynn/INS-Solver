inline void ls_flux(const int *edgeNum, const double **x, const double **y,
                    const double **sJ, const double **nx, const double **ny,
                    const double **s, double **dsldx, double **dsrdx,
                    double **dsldy, double **dsrdy) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse;

  if(edgeR == 0) {
    if(edgeL == 0) {
      reverse = !(x[0][0] == x[1][0] && y[0][0] == y[1][0]);
    } else if(edgeL == 1) {
      reverse = !(x[0][1] == x[1][0] && y[0][1] == y[1][0]);
    } else {
      reverse = !(x[0][2] == x[1][0] && y[0][2] == y[1][0]);
    }
  } else if(edgeR == 1) {
    if(edgeL == 0) {
      reverse = !(x[0][0] == x[1][1] && y[0][0] == y[1][1]);
    } else if(edgeL == 1) {
      reverse = !(x[0][1] == x[1][1] && y[0][1] == y[1][1]);
    } else {
      reverse = !(x[0][2] == x[1][1] && y[0][2] == y[1][1]);
    }
  } else {
    if(edgeL == 0) {
      reverse = !(x[0][0] == x[1][2] && y[0][0] == y[1][2]);
    } else if(edgeL == 1) {
      reverse = !(x[0][1] == x[1][2] && y[0][1] == y[1][2]);
    } else {
      reverse = !(x[0][2] == x[1][2] && y[0][2] == y[1][2]);
    }
  }

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
      rInd = exIndR + 7 + i - 1;
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

  for(int i = 0; i < 7; i++) {
    int lInd;
    int rInd = exIndR + i;
    if(reverse) {
      lInd = exIndL + 7 + i - 1;
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
