inline void poisson_coarse_matrix_2d_op1(const DG_FP *dr, const DG_FP *ds, 
                                         const DG_FP *mass, const DG_FP *rx, 
                                         const DG_FP *sx, const DG_FP *ry, 
                                         const DG_FP *sy, const DG_FP *J, 
                                         DG_FP *op1) {
  const DG_FP *dr_mat = dr;
  const DG_FP *ds_mat = ds;
  const DG_FP *mass_mat = mass;

  DG_FP Dx[DG_NP_N1 * DG_NP_N1], Dy[DG_NP_N1 * DG_NP_N1];
  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);
      Dx[ind] = rx[0] * dr_mat[ind] + sx[0] * ds_mat[ind];
      Dy[ind] = ry[0] * dr_mat[ind] + sy[0] * ds_mat[ind];
    }
  }

  DG_FP Dx_t[DG_NP_N1 * DG_NP_N1], Dy_t[DG_NP_N1 * DG_NP_N1];
  op2_in_kernel_gemm(false, false, DG_NP_N1, DG_NP_N1, DG_NP_N1, 1.0, mass_mat, DG_NP_N1, Dx, DG_NP_N1, 0.0, Dx_t, DG_NP_N1);
  op2_in_kernel_gemm(false, false, DG_NP_N1, DG_NP_N1, DG_NP_N1, 1.0, mass_mat, DG_NP_N1, Dy, DG_NP_N1, 0.0, Dy_t, DG_NP_N1);

  op2_in_kernel_gemm(true, false, DG_NP_N1, DG_NP_N1, DG_NP_N1, J[0], Dx, DG_NP_N1, Dx_t, DG_NP_N1, 0.0, op1, DG_NP_N1);
  op2_in_kernel_gemm(true, false, DG_NP_N1, DG_NP_N1, DG_NP_N1, J[0], Dy, DG_NP_N1, Dy_t, DG_NP_N1, 1.0, op1, DG_NP_N1);
}
