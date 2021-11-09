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

  for(int j = 0; j < 10; j++) {
    for(int i = 0; i < 6; i++) {
      int indL, indR;
      if(!reverse) {
        indL = j * 6 + i;
        indR = j * 6 + i;
      } else {
        indL = j * 6 + i;
        indR = j * 6 + (10 - 1 - i);
      }

      gVPL[indL] = gFL[indR];
      gVPR[indR] = gFR[indL];
    }
  }
}
