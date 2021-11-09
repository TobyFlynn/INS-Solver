inline void init_gauss_grad5(const int *edgeNum,
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
    DxR = Dx0L;
    DyR = Dy0L;
  } else if(edgeL == 1) {
    DxR = Dx1L;
    DyR = Dy1L;
  } else {
    DxR = Dx2L;
    DyR = Dy2L;
  }

  if(edgeR == 0) {
    DxL = Dx0R;
    DyL = Dy0R;
  } else if(edgeR == 1) {
    DxL = Dx1R;
    DyL = Dy1R;
  } else {
    DxL = Dx2R;
    DyL = Dy2R;
  }

  for(int m = 0; m < DG_GF_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind  = m * DG_NP + n;
      int indL = edgeL * DG_GF_NP + m;
      int indR = edgeR * DG_GF_NP + m;

      dL[ind] = nxL[indL] * DxL[ind] + nyL[indL] * DyL[ind];
      dR[ind] = nxR[indR] * DxR[ind] + nyR[indR] * DyR[ind];
    }
  }
}