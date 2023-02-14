inline void transpose_edges(DG_FP *opL, DG_FP *opR) {
  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < m; n++) {
      DG_FP tmpL = opL[m * DG_NP + n];
      opL[m * DG_NP + n] = opL[n * DG_NP + m];
      opL[n * DG_NP + m] = tmpL;

      DG_FP tmpR = opR[m * DG_NP + n];
      opR[m * DG_NP + n] = opR[n * DG_NP + m];
      opR[n * DG_NP + m] = tmpR;
    }
  }
}