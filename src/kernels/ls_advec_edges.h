inline void ls_advec_edges(const int *edgeNum, const bool *rev,
                           const double **q, double **exQ) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  // Copy data from R to L
  int exInd = edgeL * DG_NPF;
  int *fmask = &FMASK[edgeR * DG_NPF];

  for(int i = 0; i < DG_NPF; i++) {
    int rInd;
    if(reverse) {
      rInd = fmask[DG_NPF - i - 1];
    } else {
      rInd = fmask[i];
    }
    exQ[0][exInd + i] += q[1][rInd];
  }

  // Copy data from L to R
  exInd = edgeR * DG_NPF;
  fmask = &FMASK[edgeL * DG_NPF];

  for(int i = 0; i < DG_NPF; i++) {
    int lInd;
    if(reverse) {
      lInd = fmask[DG_NPF - i - 1];
    } else {
      lInd = fmask[i];
    }
    exQ[1][exInd + i] += q[0][lInd];
  }
}
