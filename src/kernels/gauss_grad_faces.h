inline void gauss_grad_faces(const int *edgeNum, const double **mDx0,
                             const double **mDy0, const double **mDx1,
                             const double **mDy1, const double **mDx2,
                             const double **mDy2, double **pDx0, double **pDy0,
                             double **pDx1, double **pDy1, double **pDx2,
                             double **pDy2) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  const double **mDxL, **mDyL, **mDxR, **mDyR;
  double **pDxL, **pDyL, **pDxR, **pDyR;
  if(edgeL == 0) {
    mDxR = mDx0;
    mDyR = mDy0;
    pDxL = pDx0;
    pDyL = pDy0;
  } else if(edgeL == 1) {
    mDxR = mDx1;
    mDyR = mDy1;
    pDxL = pDx1;
    pDyL = pDy1;
  } else {
    mDxR = mDx2;
    mDyR = mDy2;
    pDxL = pDx2;
    pDyL = pDy2;
  }

  if(edgeR == 0) {
    mDxL = mDx0;
    mDyL = mDy0;
    pDxR = pDx0;
    pDyR = pDy0;
  } else if(edgeR == 1) {
    mDxL = mDx1;
    mDyL = mDy1;
    pDxR = pDx1;
    pDyR = pDy1;
  } else {
    mDxL = mDx2;
    mDyL = mDy2;
    pDxR = pDx2;
    pDyR = pDy2;
  }

  for(int m = 0; m < DG_GF_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int indL = m * DG_NP + n;
      int indR = m * DG_NP + n;

      pDxL[0][indL] += mDxL[1][indR];
      pDyL[0][indL] += mDyL[1][indR];
      pDxR[1][indR] += mDxR[0][indL];
      pDyR[1][indR] += mDyR[0][indL];
    }
  }
}
