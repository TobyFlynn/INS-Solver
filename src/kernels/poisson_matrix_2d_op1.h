inline void poisson_matrix_2d_op1(const int *order, const DG_FP *dr,
                                  const DG_FP *ds, const DG_FP *mass,
                                  const DG_FP *rx, const DG_FP *sx,
                                  const DG_FP *ry, const DG_FP *sy,
                                  const DG_FP *J, DG_FP *op1) {
  const DG_FP *dr_mat = &dr[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *mass_mat = &mass[(*order - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];

  DG_FP Dx[DG_NP * DG_NP], Dy[DG_NP * DG_NP];
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      Dx[ind] = rx[0] * dr_mat[ind] + sx[0] * ds_mat[ind];
      Dy[ind] = ry[0] * dr_mat[ind] + sy[0] * ds_mat[ind];
    }
  }

  DG_FP Dx_t[DG_NP * DG_NP], Dy_t[DG_NP * DG_NP];
  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, 1.0, mass_mat, dg_np, Dx, dg_np, 0.0, Dx_t, dg_np);
  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, 1.0, mass_mat, dg_np, Dy, dg_np, 0.0, Dy_t, dg_np);

  op2_in_kernel_gemm(true, false, dg_np, dg_np, dg_np, J[0], Dx, dg_np, Dx_t, dg_np, 0.0, op1, dg_np);
  op2_in_kernel_gemm(true, false, dg_np, dg_np, dg_np, J[0], Dy, dg_np, Dy_t, dg_np, 1.0, op1, dg_np);
}
