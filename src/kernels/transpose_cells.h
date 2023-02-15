inline void transpose_cells(DG_FP *op) {
  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < m; n++) {
      // DG_FP tmp = op[m * DG_NP + n];
      // op[m * DG_NP + n] = op[n * DG_NP + m];
      // op[n * DG_NP + m] = tmp;

      int ind  = DG_MAT_IND(m, n, DG_NP, DG_NP);
      int indT = DG_MAT_IND(n, m, DG_NP, DG_NP);
      DG_FP tmp = op[indT];
      op[indT]  = op[ind];
      op[ind]   = tmp;
    }
  }
}
