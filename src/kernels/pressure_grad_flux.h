inline void pressure_grad_flux(const int *edgeNum, const bool *rev, const double **nx,
                               const double **ny, const double **fscale, const double **p,
                               double **pX, double **pY) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  // Copy data from R to L
  int exInd = 0;
  if(edgeL == 1) exInd = 5;
  else if(edgeL == 2) exInd = 2 * 5;

  int *fmaskL;
  if(edgeL == 0) {
    fmaskL = FMASK;
  } else if(edgeL == 1) {
    fmaskL = &FMASK[5];
  } else {
    fmaskL = &FMASK[2 * 5];
  }

  int *fmaskR;
  if(edgeR == 0) {
    fmaskR = FMASK;
  } else if(edgeR == 1) {
    fmaskR = &FMASK[5];
  } else {
    fmaskR = &FMASK[2 * 5];
  }

  for(int i = 0; i < 5; i++) {
    int lInd = fmaskL[i];
    int rInd;
    if(reverse) {
      rInd = fmaskR[5 - i - 1];
    } else {
      rInd = fmaskR[i];
    }
    double flux = p[0][lInd] - 0.5 * (p[0][lInd] + p[1][rInd]);
    pX[0][exInd + i] += fscale[0][exInd + i] * nx[0][exInd + i] * flux;
    pY[0][exInd + i] += fscale[0][exInd + i] * ny[0][exInd + i] * flux;
  }

  // Copy data from L to R
  exInd = 0;
  if(edgeR == 1) exInd = 5;
  else if(edgeR == 2) exInd = 2 * 5;

  for(int i = 0; i < 5; i++) {
    int rInd = fmaskR[i];
    int lInd;
    if(reverse) {
      lInd = fmaskL[5 - i - 1];
    } else {
      lInd = fmaskL[i];
    }
    double flux = p[1][rInd] - 0.5 * (p[0][lInd] + p[1][rInd]);
    pX[1][exInd + i] += fscale[1][exInd + i] * nx[1][exInd + i] * flux;
    pY[1][exInd + i] += fscale[1][exInd + i] * ny[1][exInd + i] * flux;
  }
}
