inline void fvmf_3d_mm_diag(const DG_FP *factor, const DG_FP *geof, DG_FP *diag) {
  const DG_FP *mass_mat = &dg_Mass_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    int op_ind = DG_MAT_IND(i, i, DG_NP, DG_NP);
    diag[i] += factor[i] * geof[J_IND] * mass_mat[op_ind];
  }
}