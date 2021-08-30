inline void gauss_gfi_faces2(const int *edgeNum, const bool *rev,
                             double *gVPL, double *gVPR) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  const double *gFL, *gFR;
  if(edgeL == 0) {
    gFR = gFInterp0_g;
  } else if(edgeL == 1) {
    gFR = gFInterp1_g;
  } else {
    gFR = gFInterp2_g;
  }

  if(edgeR == 0) {
    gFL = gFInterp0_g;
  } else if(edgeR == 1) {
    gFL = gFInterp1_g;
  } else {
    gFL = gFInterp2_g;
  }

  for(int m = 0; m < DG_GF_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int indL, indR;
      if(!reverse) {
        indL = m * DG_NP + n;
        indR = m * DG_NP + n;
      } else {
        indL = m * DG_NP + n;
        indR = (DG_GF_NP - 1 - m) * DG_NP + n;
      }

      gVPL[indL] = gFL[indR];
      gVPR[indR] = gFR[indL];
    }
  }
}
