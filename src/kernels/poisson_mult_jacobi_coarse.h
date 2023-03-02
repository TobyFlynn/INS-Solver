inline void poisson_mult_jacobi_coarse(const DG_FP *op, DG_FP *rhs) {
  for(int i = 0; i < DG_NP_N1; i++) {
    // const int op_ind = i * dg_np + i;
    int op_ind = DG_MAT_IND(i, i, DG_NP_N1, DG_NP_N1);
    rhs[i] = rhs[i] / op[op_ind];
  }
}
