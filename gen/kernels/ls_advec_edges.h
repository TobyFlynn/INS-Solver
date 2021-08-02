inline void ls_advec_edges(const int *edgeNum, const bool *rev,
                           const double **q, double **exQ) {
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
    exQ[0][exInd + i] += q[1][rInd];
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
    exQ[1][exInd + i] += q[0][lInd];
  }
}
