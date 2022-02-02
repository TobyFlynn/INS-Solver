inline void poisson_transpose_edges(double *mat0, double *mat1) {
  double tmp[DG_NP][DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      tmp[i][j] = mat0[i * DG_NP + j];
    }
  }

  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      mat0[i * DG_NP + j] = tmp[j][i];
    }
  }

  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      tmp[i][j] = mat1[i * DG_NP + j];
    }
  }

  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      mat1[i * DG_NP + j] = tmp[j][i];
    }
  }
}
