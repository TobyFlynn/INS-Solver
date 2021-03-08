inline void gauss_gfi_faces(const int *edgeNum, const double **x,
                            const double **y, double **gf0, double **gf1,
                            double **gf2) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse;

  if(edgeR == 0) {
    if(edgeL == 0) {
      reverse = !(x[0][0] == x[1][0] && y[0][0] == y[1][0]);
    } else if(edgeL == 1) {
      reverse = !(x[0][1] == x[1][0] && y[0][1] == y[1][0]);
    } else {
      reverse = !(x[0][2] == x[1][0] && y[0][2] == y[1][0]);
    }
  } else if(edgeR == 1) {
    if(edgeL == 0) {
      reverse = !(x[0][0] == x[1][1] && y[0][0] == y[1][1]);
    } else if(edgeL == 1) {
      reverse = !(x[0][1] == x[1][1] && y[0][1] == y[1][1]);
    } else {
      reverse = !(x[0][2] == x[1][1] && y[0][2] == y[1][1]);
    }
  } else {
    if(edgeL == 0) {
      reverse = !(x[0][0] == x[1][2] && y[0][0] == y[1][2]);
    } else if(edgeL == 1) {
      reverse = !(x[0][1] == x[1][2] && y[0][1] == y[1][2]);
    } else {
      reverse = !(x[0][2] == x[1][2] && y[0][2] == y[1][2]);
    }
  }

  for(int m = 0; m < 7; m++) {
    for(int n = 0; n < 15; n++) {
      int indL, indR;
      if(reverse) {
        indL = m * 15 + n;
        indR = m * 15 + n;
      } else {
        indL = m * 15 + n;
        indR = (6 - m) * 15 + n;
      }

      if(edgeL == 0) {
        if(edgeR == 0) {
          gf0[0][indL] += gFInterp0[indR];
          gf0[1][indR] += gFInterp0[indL];
        } else if(edgeR == 1) {
          gf0[0][indL] += gFInterp1[indR];
          gf1[1][indR] += gFInterp0[indL];
        } else {
          gf0[0][indL] += gFInterp2[indR];
          gf2[1][indR] += gFInterp0[indL];
        }
      } else if(edgeL == 1) {
        if(edgeR == 0) {
          gf1[0][indL] += gFInterp0[indR];
          gf0[1][indR] += gFInterp1[indL];
        } else if(edgeR == 1) {
          gf1[0][indL] += gFInterp1[indR];
          gf1[1][indR] += gFInterp1[indL];
        } else {
          gf1[0][indL] += gFInterp2[indR];
          gf2[1][indR] += gFInterp1[indL];
        }
      } else {
        if(edgeR == 0) {
          gf2[0][indL] += gFInterp0[indR];
          gf0[1][indR] += gFInterp2[indL];
        } else if(edgeR == 1) {
          gf2[0][indL] += gFInterp1[indR];
          gf1[1][indR] += gFInterp2[indL];
        } else {
          gf2[0][indL] += gFInterp2[indR];
          gf2[1][indR] += gFInterp2[indL];
        }
      }
    }
  }
}
