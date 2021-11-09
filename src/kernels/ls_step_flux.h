inline void ls_step_flux(const int *edgeNum, const bool *rev,
                         const double **fscale, const double **step,
                         double **flux) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  // Copy data from R to L
  int exInd = edgeL * DG_NPF;
  int *fmaskL = &FMASK[edgeL * DG_NPF];
  int *fmaskR = &FMASK[edgeR * DG_NPF];

  for(int i = 0; i < DG_NPF; i++) {
    int lInd = fmaskL[i];
    int rInd;
    if(reverse) {
      rInd = fmaskR[DG_NPF - i - 1];
    } else {
      rInd = fmaskR[i];
    }
    double tmp = step[0][lInd] - (step[0][lInd] + step[1][rInd]) / 2.0;
    flux[0][exInd + i] += fscale[0][exInd + i] * tmp;
  }

  // Copy data from L to R
  exInd = edgeR * DG_NPF;

  for(int i = 0; i < DG_NPF; i++) {
    int rInd = fmaskR[i];
    int lInd;
    if(reverse) {
      lInd = fmaskL[DG_NPF - i - 1];
    } else {
      lInd = fmaskL[i];
    }
    double tmp = step[1][rInd] - (step[0][lInd] + step[1][rInd]) / 2.0;
    flux[1][exInd + i] += fscale[1][exInd + i] * tmp;
  }
}
