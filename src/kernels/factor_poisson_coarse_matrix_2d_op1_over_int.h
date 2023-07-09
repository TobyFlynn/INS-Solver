inline void factor_poisson_coarse_matrix_2d_op1_over_int(const DG_FP *dr, const DG_FP *ds, 
                            const DG_FP *mass, const DG_FP *rx, const DG_FP *sx,
                            const DG_FP *ry, const DG_FP *sy, const DG_FP *J,
                            const DG_FP *factor, DG_FP *op1) {
  const DG_FP *dr_mat = dr;
  const DG_FP *ds_mat = ds;
  const DG_FP *mass_mat = mass;

  DG_FP Dx[DG_NP_N1 * DG_NP_N1], Dy[DG_NP_N1 * DG_NP_N1];
  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int ind = i + j * DG_NP_N1;
      int ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);
      Dx[ind] = rx[0] * dr_mat[ind] + sx[0] * ds_mat[ind];
      Dy[ind] = ry[0] * dr_mat[ind] + sy[0] * ds_mat[ind];
    }
  }

  DG_FP Dx_t[DG_NP_N1 * DG_NP_N1], Dy_t[DG_NP_N1 * DG_NP_N1];
  op2_in_kernel_gemm(false, false, DG_NP_N1, DG_NP_N1, DG_NP_N1, 1.0, mass_mat, DG_NP_N1, Dx, DG_NP_N1, 0.0, Dx_t, DG_NP_N1);
  op2_in_kernel_gemm(false, false, DG_NP_N1, DG_NP_N1, DG_NP_N1, 1.0, mass_mat, DG_NP_N1, Dy, DG_NP_N1, 0.0, Dy_t, DG_NP_N1);

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
        DG_FP tmp = Dx[a_ind] * Dx_t[b_ind] + Dy[a_ind] * Dy_t[b_ind];
        op1[op_ind] += factor[k] * tmp;
      }
      op1[op_ind] *= J[0];
    }
  }
}
