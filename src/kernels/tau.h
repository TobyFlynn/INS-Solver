inline void tau(const int *edgeNum, const double **J, const double **sJ, double **tau) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  double hminL = 2.0 * J[0][FMASK[edgeL * 5]] / sJ[0][edgeL * 5];
  double hminR = 2.0 * J[1][FMASK[edgeR * 5]] / sJ[1][edgeR * 5];

  if(hminL > hminR) {
    tau[0][edgeL] += 100.0 * 15.0 / hminR;
    tau[1][edgeR] += 100.0 * 15.0 / hminR;
  } else {
    tau[0][edgeL] += 100.0 * 15.0 / hminL;
    tau[1][edgeR] += 100.0 * 15.0 / hminL;
  }
}
