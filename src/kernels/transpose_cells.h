inline void transpose_cells(double *op) {
  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < m; n++) {
      double tmp = op[m * DG_NP + n];
      op[m * DG_NP + n] = op[n * DG_NP + m];
      op[n * DG_NP + m] = tmp;
    }
  }
}