inline void poisson_matrix_2d_op1_diag(const int *order, const DG_FP *geof,
                                       DG_FP *diag) {
  const DG_FP *dr_mat = dg_Dr_kernel + (*order - 1) * DG_NP * DG_NP;
  const DG_FP *ds_mat = dg_Ds_kernel + (*order - 1) * DG_NP * DG_NP;
  const DG_FP *mass_mat = dg_Mass_kernel + (*order - 1) * DG_NP * DG_NP;
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];

  DG_FP Dx[DG_NP * DG_NP], Dy[DG_NP * DG_NP];
  const DG_FP rx = geof[RX_IND];
  const DG_FP ry = geof[RY_IND];
  const DG_FP sx = geof[SX_IND];
  const DG_FP sy = geof[SY_IND];
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      Dx[ind] = rx * dr_mat[ind] + sx * ds_mat[ind];
      Dy[ind] = ry * dr_mat[ind] + sy * ds_mat[ind];
    }
  }

  DG_FP Dx_t[DG_NP * DG_NP], Dy_t[DG_NP * DG_NP];
  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, 1.0, mass_mat, dg_np, Dx, dg_np, 0.0, Dx_t, dg_np);
  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, 1.0, mass_mat, dg_np, Dy, dg_np, 0.0, Dy_t, dg_np);

  const DG_FP J = geof[J_IND];

  for(int i = 0; i < dg_np; i++) {
    diag[i] = 0.0;
    for(int k = 0; k < dg_np; k++) {
      // int a_ind = i * dg_np + k;
      int a_ind = DG_MAT_IND(k, i, dg_np, dg_np);
      // int b_ind = i * dg_np + k;
      int b_ind = DG_MAT_IND(k, i, dg_np, dg_np);
      diag[i] += Dx[a_ind] * Dx_t[b_ind] + Dy[a_ind] * Dy_t[b_ind];
    }
    diag[i] *= J;
  }
}
