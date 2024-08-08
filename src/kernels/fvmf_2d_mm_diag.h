inline void fvmf_2d_mm_diag(const DG_FP *factor, const DG_FP *geof, DG_FP *diag) {
  const DG_FP *mass_mat = &dg_Mass_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];

  const DG_FP J = geof[J_IND];
  for(int i = 0; i < DG_NP; i++) {
    int op_ind = DG_MAT_IND(i, i, DG_NP, DG_NP);
    diag[i] += factor[i] * J * mass_mat[op_ind];
  }
}