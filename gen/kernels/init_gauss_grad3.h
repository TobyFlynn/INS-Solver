inline void init_gauss_grad3(const int *edgeNum,
                             const double *nxL, const double *nxR,
                             const double *nyL, const double *nyR,
                             const double *Dx0L, const double *Dx0R,
                             const double *Dy0L, const double *Dy0R,
                             const double *Dx1L, const double *Dx1R,
                             const double *Dy1L, const double *Dy1R,
                             const double *Dx2L, const double *Dx2R,
                             const double *Dy2L, const double *Dy2R,
                             double *dL, double *dR) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  const double *DxL, *DyL, *DxR, *DyR;

  if(edgeL == 0) {
    DxL = Dx0L;
    DyL = Dy0L;
  } else if(edgeL == 1) {
    DxL = Dx1L;
    DyL = Dy1L;
  } else {
    DxL = Dx2L;
    DyL = Dy2L;
  }

  if(edgeR == 0) {
    DxR = Dx0R;
    DyR = Dy0R;
  } else if(edgeR == 1) {
    DxR = Dx1R;
    DyR = Dy1R;
  } else {
    DxR = Dx2R;
    DyR = Dy2R;
  }

  for(int m = 0; m < 6; m++) {
    for(int n = 0; n < 10; n++) {
      int ind  = m * 10 + n;
      int indL = edgeL * 6 + m;
      int indR = edgeR * 6 + m;

      dL[ind] = nxL[indL] * DxL[ind] + nyL[indL] * DyL[ind];
      dR[ind] = nxR[indR] * DxR[ind] + nyR[indR] * DyR[ind];
    }
  }
}
