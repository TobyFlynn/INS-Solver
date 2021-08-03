inline void gauss_grad_faces(const int *edgeNum, const double **mDx0,
                             const double **mDy0, const double **mDx1,
                             const double **mDy1, const double **mDx2,
                             const double **mDy2, double **pDx0, double **pDy0,
                             double **pDx1, double **pDy1, double **pDx2,
                             double **pDy2) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  for(int m = 0; m < 6; m++) {
    for(int n = 0; n < 10; n++) {
      int indL = m * 10 + n;
      int indR = m * 10 + n;

      if(edgeL == 0) {
        if(edgeR == 0) {
          pDx0[0][indL] += mDx0[1][indR];
          pDy0[0][indL] += mDy0[1][indR];
          pDx0[1][indR] += mDx0[0][indL];
          pDy0[1][indR] += mDy0[0][indL];
        } else if(edgeR == 1) {
          pDx0[0][indL] += mDx1[1][indR];
          pDy0[0][indL] += mDy1[1][indR];
          pDx1[1][indR] += mDx0[0][indL];
          pDy1[1][indR] += mDy0[0][indL];
        } else {
          pDx0[0][indL] += mDx2[1][indR];
          pDy0[0][indL] += mDy2[1][indR];
          pDx2[1][indR] += mDx0[0][indL];
          pDy2[1][indR] += mDy0[0][indL];
        }
      } else if(edgeL == 1) {
        if(edgeR == 0) {
          pDx1[0][indL] += mDx0[1][indR];
          pDy1[0][indL] += mDy0[1][indR];
          pDx0[1][indR] += mDx1[0][indL];
          pDy0[1][indR] += mDy1[0][indL];
        } else if(edgeR == 1) {
          pDx1[0][indL] += mDx1[1][indR];
          pDy1[0][indL] += mDy1[1][indR];
          pDx1[1][indR] += mDx1[0][indL];
          pDy1[1][indR] += mDy1[0][indL];
        } else {
          pDx1[0][indL] += mDx2[1][indR];
          pDy1[0][indL] += mDy2[1][indR];
          pDx2[1][indR] += mDx1[0][indL];
          pDy2[1][indR] += mDy1[0][indL];
        }
      } else {
        if(edgeR == 0) {
          pDx2[0][indL] += mDx0[1][indR];
          pDy2[0][indL] += mDy0[1][indR];
          pDx0[1][indR] += mDx2[0][indL];
          pDy0[1][indR] += mDy2[0][indL];
        } else if(edgeR == 1) {
          pDx2[0][indL] += mDx1[1][indR];
          pDy2[0][indL] += mDy1[1][indR];
          pDx1[1][indR] += mDx2[0][indL];
          pDy1[1][indR] += mDy2[0][indL];
        } else {
          pDx2[0][indL] += mDx2[1][indR];
          pDy2[0][indL] += mDy2[1][indR];
          pDx2[1][indR] += mDx2[0][indL];
          pDy2[1][indR] += mDy2[0][indL];
        }
      }
    }
  }
}
