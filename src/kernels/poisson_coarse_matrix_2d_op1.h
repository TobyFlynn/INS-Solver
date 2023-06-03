inline void poisson_coarse_matrix_2d_op1(const DG_FP *geof, DG_FP *op1) {
  const DG_FP *dr_mat = dg_Dr_kernel;
  const DG_FP *ds_mat = dg_Ds_kernel;
  const DG_FP *mass_mat = dg_Mass_kernel;

  DG_FP Dx[DG_NP_N1 * DG_NP_N1], Dy[DG_NP_N1 * DG_NP_N1];
  const DG_FP rx = geof[RX_IND];
  const DG_FP ry = geof[RY_IND];
  const DG_FP sx = geof[SX_IND];
  const DG_FP sy = geof[SY_IND];
  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);
      Dx[ind] = rx * dr_mat[ind] + sx * ds_mat[ind];
      Dy[ind] = ry * dr_mat[ind] + sy * ds_mat[ind];
    }
  }

  DG_FP Dx_t[DG_NP_N1 * DG_NP_N1], Dy_t[DG_NP_N1 * DG_NP_N1];
  op2_in_kernel_gemm(false, false, DG_NP_N1, DG_NP_N1, DG_NP_N1, 1.0, mass_mat, DG_NP_N1, Dx, DG_NP_N1, 0.0, Dx_t, DG_NP_N1);
  op2_in_kernel_gemm(false, false, DG_NP_N1, DG_NP_N1, DG_NP_N1, 1.0, mass_mat, DG_NP_N1, Dy, DG_NP_N1, 0.0, Dy_t, DG_NP_N1);

  const DG_FP J = geof[J_IND];
  op2_in_kernel_gemm(true, false, DG_NP_N1, DG_NP_N1, DG_NP_N1, J, Dx, DG_NP_N1, Dx_t, DG_NP_N1, 0.0, op1, DG_NP_N1);
  op2_in_kernel_gemm(true, false, DG_NP_N1, DG_NP_N1, DG_NP_N1, J, Dy, DG_NP_N1, Dy_t, DG_NP_N1, 1.0, op1, DG_NP_N1);
}
