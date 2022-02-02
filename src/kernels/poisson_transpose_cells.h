inline void poisson_transpose_cells(double *mat) {
  double tmp[DG_NP][DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      tmp[i][j] = mat[i * DG_NP + j];
    }
  }

  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      mat[i * DG_NP + j] = tmp[j][i];
    }
  }
}
