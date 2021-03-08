inline void gauss_grad_faces(const int *edgeNum, const double **x,
                             const double **y, const double **mDx0,
                             const double **mDy0, const double **mDx1,
                             const double **mDy1, const double **mDx2,
                             const double **mDy2, double **pDx0, double **pDy0,
                             double **pDx1, double **pDy1, double **pDx2,
                             double **pDy2) {
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

  for(int m = 0; m < 7; m++) {
    for(int n = 0; n < 15; n++) {
      int indL, indR;
      if(reverse) {
        indL = m * 15 + n;
        indR = m * 15 + n;
      } else {
        indL = m * 15 + n;
        indR = (6 - m) * 15 + n;
      }

      if(edgeL == 0) {
        if(edgeR == 0) {
          pDx0[0][indL] += mDx0[1][indR];
          pDy0[0][indL] += mDy0[1][indR];
          pDx0[1][indR] += mDx0[0][indL];
          pDx0[1][indR] += mDx0[0][indL];
        } else if(edgeR == 1) {
          pDx0[0][indL] += mDx1[1][indR];
          pDy0[0][indL] += mDy1[1][indR];
          pDx1[1][indR] += mDx0[0][indL];
          pDx1[1][indR] += mDx0[0][indL];
        } else {
          pDx0[0][indL] += mDx2[1][indR];
          pDy0[0][indL] += mDy2[1][indR];
          pDx2[1][indR] += mDx0[0][indL];
          pDx2[1][indR] += mDx0[0][indL];
        }
      } else if(edgeL == 1) {
        if(edgeR == 0) {
          pDx1[0][indL] += mDx0[1][indR];
          pDy1[0][indL] += mDy0[1][indR];
          pDx0[1][indR] += mDx1[0][indL];
          pDx0[1][indR] += mDx1[0][indL];
        } else if(edgeR == 1) {
          pDx1[0][indL] += mDx1[1][indR];
          pDy1[0][indL] += mDy1[1][indR];
          pDx1[1][indR] += mDx1[0][indL];
          pDx1[1][indR] += mDx1[0][indL];
        } else {
          pDx1[0][indL] += mDx2[1][indR];
          pDy1[0][indL] += mDy2[1][indR];
          pDx2[1][indR] += mDx1[0][indL];
          pDx2[1][indR] += mDx1[0][indL];
        }
      } else {
        if(edgeR == 0) {
          pDx2[0][indL] += mDx0[1][indR];
          pDy2[0][indL] += mDy0[1][indR];
          pDx0[1][indR] += mDx2[0][indL];
          pDx0[1][indR] += mDx2[0][indL];
        } else if(edgeR == 1) {
          pDx2[0][indL] += mDx1[1][indR];
          pDy2[0][indL] += mDy1[1][indR];
          pDx1[1][indR] += mDx2[0][indL];
          pDx1[1][indR] += mDx2[0][indL];
        } else {
          pDx2[0][indL] += mDx2[1][indR];
          pDy2[0][indL] += mDy2[1][indR];
          pDx2[1][indR] += mDx2[0][indL];
          pDx2[1][indR] += mDx2[0][indL];
        }
      }
    }
  }
}
