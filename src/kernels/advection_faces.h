inline void advection_faces(const int *edgeNum, const bool *rev, const double **q0,
                            const double **q1, double **exQ0, double **exQ1) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  // Copy data from R to L
  int exInd = 0;
  if(edgeL == 1) exInd = DG_NPF;
  else if(edgeL == 2) exInd = 2 * DG_NPF;

  int *fmask;

  if(edgeR == 0) {
    fmask = FMASK;
  } else if(edgeR == 1) {
    fmask = &FMASK[DG_NPF];
  } else {
    fmask = &FMASK[2 * DG_NPF];
  }

  for(int i = 0; i < DG_NPF; i++) {
    int rInd;
    if(reverse) {
      rInd = fmask[DG_NPF - i - 1];
    } else {
      rInd = fmask[i];
    }
    exQ0[0][exInd + i] += q0[1][rInd];
    exQ1[0][exInd + i] += q1[1][rInd];
  }

  // Copy data from L to R
  exInd = 0;
  if(edgeR == 1) exInd = DG_NPF;
  else if(edgeR == 2) exInd = 2 * DG_NPF;

  if(edgeL == 0) {
    fmask = FMASK;
  } else if(edgeL == 1) {
    fmask = &FMASK[DG_NPF];
  } else {
    fmask = &FMASK[2 * DG_NPF];
  }

  for(int i = 0; i < DG_NPF; i++) {
    int lInd;
    if(reverse) {
      lInd = fmask[DG_NPF - i - 1];
    } else {
      lInd = fmask[i];
    }
    exQ0[1][exInd + i] += q0[0][lInd];
    exQ1[1][exInd + i] += q1[0][lInd];
  }
}
