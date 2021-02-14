inline void set_tau(const int *edgeNum, const double **x, const double **y,
                    const double **J, const double **sJ, double **tau) {
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

  // Copy data from R to L
  int exIndL = 0;
  if(edgeL == 1) exIndL = 5;
  else if(edgeL == 2) exIndL = 2 * 5;

  int exIndR = 0;
  if(edgeR == 1) exIndR = 5;
  else if(edgeR == 2) exIndR = 2 * 5;

  int *fmaskR;

  if(edgeR == 0) {
    fmaskR = FMASK;
  } else if(edgeR == 1) {
    fmaskR = &FMASK[5];
  } else {
    fmaskR = &FMASK[2 * 5];
  }

  int *fmaskL;

  if(edgeL == 0) {
    fmaskL = FMASK;
  } else if(edgeL == 1) {
    fmaskL = &FMASK[5];
  } else {
    fmaskL = &FMASK[2 * 5];
  }

  for(int i = 0; i < 5; i++) {
    int rIndF, lIndF, rInd, lInd;
    if(reverse) {
      rIndF = fmaskR[5 - i - 1];
      rInd = exIndR + 5 - i - 1;
    } else {
      rIndF = fmaskR[i];
      rInd = exIndR + i;
    }
    lIndF = fmaskL[i];
    lInd = exIndL + i;

    double lH = 2.0 * J[0][lIndF] / sJ[0][lInd];
    double rH = 2.0 * J[1][rIndF] / sJ[1][rInd];
    if(lH < rH) {
      tau[0][lInd] += 15.0 / lH;
      tau[1][rInd] += 15.0 / lH;
    } else {
      tau[0][lInd] += 15.0 / rH;
      tau[1][rInd] += 15.0 / rH;
    }
  }
}
