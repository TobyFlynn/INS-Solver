inline void ls_advec_edges(const int **p, const int *edgeNum, const bool *rev,
                           const double **q, double **exQ) {
  // Work out which edge for each element
  const int edgeL = edgeNum[0];
  const int edgeR = edgeNum[1];
  const bool reverse = *rev;

  // Get constants for these elements' orders
  const int dg_npfL  = DG_CONSTANTS[(p[0][0] - 1) * 5 + 1];
  const int dg_npfR  = DG_CONSTANTS[(p[1][0] - 1) * 5 + 1];
  const int *fmaskL_ = &FMASK[(p[0][0] - 1) * 3 * DG_NPF];
  const int *fmaskR_ = &FMASK[(p[1][0] - 1) * 3 * DG_NPF];

  // Copy data from R to L
  const int exIndL = edgeL * dg_npfL;
  const int *fmask = &fmaskR_[edgeR * dg_npfR];

  for(int i = 0; i < dg_npfR; i++) {
    int rInd;
    if(reverse) {
      rInd = fmask[dg_npfR - i - 1];
    } else {
      rInd = fmask[i];
    }
    exQ[0][exIndL + i] += q[1][rInd];
  }

  // Copy data from L to R
  const int exIndR = edgeR * dg_npfR;
  fmask = &fmaskL_[edgeL * dg_npfL];

  for(int i = 0; i < dg_npfL; i++) {
    int lInd;
    if(reverse) {
      lInd = fmask[dg_npfL - i - 1];
    } else {
      lInd = fmask[i];
    }
    exQ[1][exIndR + i] += q[0][lInd];
  }
}
