inline void factor_poisson_matrix_2d_mm_diag(const int *order,
                          const DG_FP *factor, const DG_FP *geof, DG_FP *diag) {
  const DG_FP *mass_mat = &dg_Mass_kernel[(*order - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];

  const DG_FP J = geof[J_IND];
  for(int i = 0; i < dg_np; i++) {
    int op_ind = DG_MAT_IND(i, i, dg_np, dg_np);
    diag[i] += factor[i] * J * mass_mat[op_ind];
  }
}
