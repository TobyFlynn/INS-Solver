inline void gauss_gfi_faces(const int *edgeNum, const bool *rev,
                            double **gf0, double **gf1, double **gf2) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  for(int m = 0; m < 7; m++) {
    for(int n = 0; n < 15; n++) {
      int indL, indR;
      if(!reverse) {
        indL = m * 15 + n;
        indR = m * 15 + n;
      } else {
        indL = m * 15 + n;
        indR = (7 - 1 - m) * 15 + n;
      }

      if(edgeL == 0) {
        if(edgeR == 0) {
          gf0[0][indL] += gFInterp0_g[indR];
          gf0[1][indR] += gFInterp0_g[indL];
        } else if(edgeR == 1) {
          gf0[0][indL] += gFInterp1_g[indR];
          gf1[1][indR] += gFInterp0_g[indL];
        } else {
          gf0[0][indL] += gFInterp2_g[indR];
          gf2[1][indR] += gFInterp0_g[indL];
        }
      } else if(edgeL == 1) {
        if(edgeR == 0) {
          gf1[0][indL] += gFInterp0_g[indR];
          gf0[1][indR] += gFInterp1_g[indL];
        } else if(edgeR == 1) {
          gf1[0][indL] += gFInterp1_g[indR];
          gf1[1][indR] += gFInterp1_g[indL];
        } else {
          gf1[0][indL] += gFInterp2_g[indR];
          gf2[1][indR] += gFInterp1_g[indL];
        }
      } else {
        if(edgeR == 0) {
          gf2[0][indL] += gFInterp0_g[indR];
          gf0[1][indR] += gFInterp2_g[indL];
        } else if(edgeR == 1) {
          gf2[0][indL] += gFInterp1_g[indR];
          gf1[1][indR] += gFInterp2_g[indL];
        } else {
          gf2[0][indL] += gFInterp2_g[indR];
          gf2[1][indR] += gFInterp2_g[indL];
        }
      }
    }
  }
}
