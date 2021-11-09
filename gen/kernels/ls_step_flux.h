inline void ls_step_flux(const int *edgeNum, const bool *rev,
                         const double **fscale, const double **step,
                         double **flux) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  // Copy data from R to L
  int exInd = edgeL * 4;
  int *fmaskL = &FMASK[edgeL * 4];
  int *fmaskR = &FMASK[edgeR * 4];

  for(int i = 0; i < 4; i++) {
    int lInd = fmaskL[i];
    int rInd;
    if(reverse) {
      rInd = fmaskR[4 - i - 1];
    } else {
      rInd = fmaskR[i];
    }
    double tmp = step[0][lInd] - (step[0][lInd] + step[1][rInd]) / 2.0;
    flux[0][exInd + i] += fscale[0][exInd + i] * tmp;
  }

  // Copy data from L to R
  exInd = edgeR * 4;

  for(int i = 0; i < 4; i++) {
    int rInd = fmaskR[i];
    int lInd;
    if(reverse) {
      lInd = fmaskL[4 - i - 1];
    } else {
      lInd = fmaskL[i];
    }
    double tmp = step[1][rInd] - (step[0][lInd] + step[1][rInd]) / 2.0;
    flux[1][exInd + i] += fscale[1][exInd + i] * tmp;
  }
}
