inline void pressure_grad_flux(const int *edgeNum, const bool *rev, const double **nx,
                               const double **ny, const double **fscale, const double **p,
                               double **pX, double **pY) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  // Copy data from R to L
  int exInd = edgeL * 3;
  int *fmaskL = &FMASK[edgeL * 3];
  int *fmaskR = &FMASK[edgeR * 3];

  for(int i = 0; i < 3; i++) {
    int lInd = fmaskL[i];
    int rInd;
    if(reverse) {
      rInd = fmaskR[3 - i - 1];
    } else {
      rInd = fmaskR[i];
    }
    double flux = p[0][lInd] - 0.5 * (p[0][lInd] + p[1][rInd]);
    pX[0][exInd + i] += fscale[0][exInd + i] * nx[0][exInd + i] * flux;
    pY[0][exInd + i] += fscale[0][exInd + i] * ny[0][exInd + i] * flux;
  }

  // Copy data from L to R
  exInd = edgeR * 3;

  for(int i = 0; i < 3; i++) {
    int rInd = fmaskR[i];
    int lInd;
    if(reverse) {
      lInd = fmaskL[3 - i - 1];
    } else {
      lInd = fmaskL[i];
    }
    double flux = p[1][rInd] - 0.5 * (p[0][lInd] + p[1][rInd]);
    pX[1][exInd + i] += fscale[1][exInd + i] * nx[1][exInd + i] * flux;
    pY[1][exInd + i] += fscale[1][exInd + i] * ny[1][exInd + i] * flux;
  }
}
