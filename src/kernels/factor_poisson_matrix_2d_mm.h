inline void factor_poisson_matrix_2d_mm(const int *order, const DG_FP *factor,
                                        const DG_FP *mass, const DG_FP *J,
                                        DG_FP *op1) {
  const DG_FP *mass_mat = &mass[(*order - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int op_ind = i + j * dg_np;
      int op_ind = DG_MAT_IND(i, j, dg_np, dg_np);
      op1[op_ind] += factor[i] * J[0] * mass_mat[op_ind];
    }
  }
}
