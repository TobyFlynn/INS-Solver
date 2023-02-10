inline void poisson_matrix_3d_apply_bc(const double *op, const double *bc,
                                       double *rhs) {
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NPF; j++) {
      int ind = i + j * DG_NP;
      rhs[i] += op[ind] * bc[j];
    }
  }
}
