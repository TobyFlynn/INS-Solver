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

  for(int i = 0; i < 6 * 10; i++) {
    gVPL[i] = 0.0;
    gVPR[i] = 0.0;
  }

  for(int m = 0; m < 6; m++) {
    for(int n = 0; n < 10; n++) {
      int indL, indR;
      if(!reverse) {
        indL = m * 10 + n;
        indR = m * 10 + n;
      } else {
        indL = m * 10 + n;
        indR = (6 - 1 - m) * 10 + n;
      }

      gVPL[indL] += gFL[indR];
      gVPR[indR] += gFR[indL];
    }
  }
}
