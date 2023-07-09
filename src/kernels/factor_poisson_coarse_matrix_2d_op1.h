inline void factor_poisson_coarse_matrix_2d_op1(const DG_FP *geof,
                                            const DG_FP *factor, DG_FP *op1) {
  const DG_FP *dr_mat = dg_Dr_kernel;
  const DG_FP *ds_mat = dg_Ds_kernel;
  const DG_FP *mass_mat = dg_Mass_kernel;

  DG_FP Dx[DG_NP_N1 * DG_NP_N1], Dy[DG_NP_N1 * DG_NP_N1];
  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int ind = i + j * DG_NP_N1;
      int ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);
      Dx[ind] = geof[RX_IND] * dr_mat[ind] + geof[SX_IND] * ds_mat[ind];
      Dy[ind] = geof[RY_IND] * dr_mat[ind] + geof[SY_IND] * ds_mat[ind];
    }
  }

  DG_FP Dx_t[DG_NP_N1 * DG_NP_N1], Dy_t[DG_NP_N1 * DG_NP_N1];
  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int op_ind = i + j * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);
      Dx_t[op_ind] = 0.0;
      Dy_t[op_ind] = 0.0;
      for(int k = 0; k < DG_NP_N1; k++) {
        // int a_ind = i + k * DG_NP_N1;
        int a_ind = DG_MAT_IND(i, k, DG_NP_N1, DG_NP_N1);
        // int b_ind = j * DG_NP_N1 + k;
        int b_ind = DG_MAT_IND(k, j, DG_NP_N1, DG_NP_N1);
        Dx_t[op_ind] += mass_mat[a_ind] * factor[k] * Dx[b_ind];
        Dy_t[op_ind] += mass_mat[a_ind] * factor[k] * Dy[b_ind];
      }
    }
  }

  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int op_ind = i + j * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);
      op1[op_ind] = 0.0;
      for(int k = 0; k < DG_NP_N1; k++) {
        // int a_ind = i * DG_NP_N1 + k;
        int a_ind = DG_MAT_IND(k, i, DG_NP_N1, DG_NP_N1);
        // int b_ind = j * DG_NP_N1 + k;
        int b_ind = DG_MAT_IND(k, j, DG_NP_N1, DG_NP_N1);
        op1[op_ind] += Dx[a_ind] * Dx_t[b_ind] + Dy[a_ind] * Dy_t[b_ind];
      }
      op1[op_ind] *= geof[J_IND];
    }
  }
}
