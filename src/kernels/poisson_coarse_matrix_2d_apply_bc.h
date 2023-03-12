inline void poisson_coarse_matrix_2d_apply_bc(const int *bedgeNum, const DG_FP *op, 
                                              const DG_FP *bc, DG_FP *rhs) {
  int exInd = *bedgeNum * DG_GF_NP;
  for(int m = 0; m < DG_NP_N1; m++) {
    for(int n = 0; n < DG_GF_NP; n++) {
      // int ind = m + n * dg_np;
      int ind = DG_MAT_IND(m, n, DG_NP_N1, DG_GF_NP);
      rhs[m] += op[ind] * bc[exInd + n];
    }
  }
}
