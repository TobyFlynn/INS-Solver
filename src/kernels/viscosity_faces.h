inline void viscosity_faces(const int *edgeNum, const double **q0, const double **q1, double **exQ0, double **exQ1) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  // Copy data from R to L
  int exInd = 0;
  if(edgeL == 1) exInd = 5;
  else if(edgeL == 2) exInd = 2 * 5;

  int *fmask;

  if(edgeL == 0) {
    fmask = FMASK;
  } else if(edgeL == 1) {
    fmask = &FMASK[5];
  } else {
    fmask = &FMASK[2 * 5];
  }

  for(int i = 0; i < 5; i++) {
    int lInd = fmask[i];
    exQ0[0][exInd + i] += q0[0][lInd];
    exQ1[0][exInd + i] += q1[0][lInd];
  }

  // Copy data from L to R
  exInd = 0;
  if(edgeR == 1) exInd = 5;
  else if(edgeR == 2) exInd = 2 * 5;

  if(edgeR == 0) {
    fmask = FMASK;
  } else if(edgeR == 1) {
    fmask = &FMASK[5];
  } else {
    fmask = &FMASK[2 * 5];
  }

  for(int i = 0; i < 5; i++) {
    int rInd = fmask[i];
    exQ0[1][exInd + i] += q0[1][rInd];
    exQ1[1][exInd + i] += q1[1][rInd];
  }
}
