inline void fvmf_3d_mm_block_diag(const DG_FP *factor, const DG_FP *geof, DG_FP *op1) {
  const DG_FP *mass_mat = &dg_Mass_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      // int op_ind = i + j * DG_NP;
      int op_ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      op1[op_ind] += factor[i] * geof[J_IND] * mass_mat[op_ind];
    }
  }
}