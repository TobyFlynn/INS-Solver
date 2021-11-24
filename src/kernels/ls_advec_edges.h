inline void ls_advec_edges(const int **p, const int *edgeNum, const bool *rev,
                           const double **q, double **exQ) {
  // Work out which edge for each element
  const int edgeL = edgeNum[0];
  const int edgeR = edgeNum[1];
  const bool reverse = *rev;

  const int exIndL = edgeL * DG_GF_NP;
  const int exIndR = edgeR * DG_GF_NP;

  // Copy data from R to L
  for(int i = 0; i < DG_GF_NP; i++) {
    int rInd;
    if(reverse) {
      rInd = exIndR + DG_GF_NP - 1 - i;
    } else {
      rInd = exIndR + i;
    }
    exQ[0][exIndL + i] += q[1][rInd];
  }

  // Copy data from L to R
  for(int i = 0; i < DG_GF_NP; i++) {
    int lInd;
    if(reverse) {
      lInd = exIndL + DG_GF_NP - 1 - i;
    } else {
      lInd = exIndL + i;
    }
    exQ[1][exIndR + i] += q[0][lInd];
  }
}
