inline void transpose_edges(double *opL, double *opR) {
  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < m; n++) {
      double tmpL = opL[m * DG_NP + n];
      opL[m * DG_NP + n] = opL[n * DG_NP + m];
      opL[n * DG_NP + m] = tmpL;

      double tmpR = opR[m * DG_NP + n];
      opR[m * DG_NP + n] = opR[n * DG_NP + m];
      opR[n * DG_NP + m] = tmpR;
    }
  }
}
