inline void gauss_gfi_faces(const int *edgeNum, const bool *rev,
                            double **gf0, double **gf1, double **gf2) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

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
