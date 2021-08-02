inline void advection_faces(const int *edgeNum, const bool *rev, const double **q0,
                            const double **q1, double **exQ0, double **exQ1) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  // Copy data from R to L
  int exInd = edgeL * 4;
  int *fmask = &FMASK[edgeR * 4];

  for(int i = 0; i < 4; i++) {
    int rInd;
    if(reverse) {
      rInd = fmask[4 - i - 1];
    } else {
      rInd = fmask[i];
    }
    exQ0[0][exInd + i] += q0[1][rInd];
    exQ1[0][exInd + i] += q1[1][rInd];
  }

  // Copy data from L to R
  exInd = edgeR * 4;
  fmask = &FMASK[edgeL * 4];

  for(int i = 0; i < 4; i++) {
    int lInd;
    if(reverse) {
      lInd = fmask[4 - i - 1];
    } else {
      lInd = fmask[i];
    }
    exQ0[1][exInd + i] += q0[0][lInd];
    exQ1[1][exInd + i] += q1[0][lInd];
  }
}
