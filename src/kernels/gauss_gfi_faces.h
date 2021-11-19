inline void gauss_gfi_faces(const int *edgeNum, const bool *rev,
                            double **gf0, double **gf1, double **gf2) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  const double *gVML, *gVMR;
  double **gFL, **gFR;
  if(edgeL == 0) {
    gVMR = gFInterp0_g;
    gFL  = gf0;
  } else if(edgeL == 1) {
    gVMR = gFInterp1_g;
    gFL  = gf1;
  } else {
    gVMR = gFInterp2_g;
    gFL  = gf2;
  }

  if(edgeR == 0) {
    gVML = gFInterp0_g;
    gFR  = gf0;
  } else if(edgeR == 1) {
    gVML = gFInterp1_g;
    gFR  = gf1;
  } else {
    gVML = gFInterp2_g;
    gFR  = gf2;
  }

  for(int m = 0; m < DG_GF_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int indL, indR;
      int indL_col, indR_col;
      if(!reverse) {
        indL = m * DG_NP + n;
        indR = m * DG_NP + n;
        indL_col = m + n * DG_GF_NP;
        indR_col = m + n * DG_GF_NP;
      } else {
        indL = m * DG_NP + n;
        indR = (DG_GF_NP - 1 - m) * DG_NP + n;
        indL_col = m + n * DG_GF_NP;
        indR_col = (DG_GF_NP - 1 - m) + n * DG_GF_NP;
      }

      gFL[0][indL] += gVML[indR_col];
      gFR[1][indR] += gVMR[indL_col];
    }
  }
}
