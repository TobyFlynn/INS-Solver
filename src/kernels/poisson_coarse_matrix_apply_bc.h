inline void poisson_coarse_matrix_apply_bc(const DG_FP *op, const DG_FP *bc,
                                           DG_FP *rhs) {
  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NPF_N1; j++) {
      // int ind = i + j * DG_NP;
      int ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NPF_N1);
      rhs[i] += op[ind] * bc[j];
    }
  }
}
