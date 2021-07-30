inline void gauss_tau(const int *edgeNum, const double **fscale, double **tau) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  if(fscale[0][edgeL * 3] > fscale[1][edgeR * 3]) {
    tau[0][edgeL] += 20 * 25 * fscale[0][edgeL * 3];
    tau[1][edgeR] += 20 * 25 * fscale[0][edgeL * 3];
  } else {
    tau[0][edgeL] += 20 * 25 * fscale[1][edgeR * 3];
    tau[1][edgeR] += 20 * 25 * fscale[1][edgeR * 3];
  }
}
