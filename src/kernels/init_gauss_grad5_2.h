inline void init_gauss_grad5_2(const int *edgeNum, const bool *rev,
                             const double *nxL, const double *nxR,
                             const double *nyL, const double *nyR,
                             const double *Dx0L, const double *Dx0R,
                             const double *Dy0L, const double *Dy0R,
                             const double *Dx1L, const double *Dx1R,
                             const double *Dy1L, const double *Dy1R,
                             const double *Dx2L, const double *Dx2R,
                             const double *Dy2L, const double *Dy2R,
                             const double *factL, const double *factR,
                             double *dL, double *dR) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

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
      int factLInd, factRInd;
      if(reverse) {
        factLInd = edgeL * DG_GF_NP + DG_GF_NP - 1 - m;
        factRInd = edgeR * DG_GF_NP + DG_GF_NP - 1 - m;
      } else {
        factLInd = edgeL * DG_GF_NP + m;
        factRInd = edgeR * DG_GF_NP + m;
      }

      dL[ind] = nxL[indL] * factR[factRInd] * DxL[ind] + nyL[indL] * factR[factRInd] * DyL[ind];
      dR[ind] = nxR[indR] * factL[factLInd] * DxR[ind] + nyR[indR] * factL[factLInd] * DyR[ind];

      // dL[ind] = nxL[indL] * factL[indL] * DxL[ind] + nyL[indL] * factL[indL] * DyL[ind];
      // dR[ind] = nxR[indR] * factR[indR] * DxR[ind] + nyR[indR] * factR[indR] * DyR[ind];
    }
  }
}
