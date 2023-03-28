inline void factor_poisson_matrix_3d_op1(const int *order, const DG_FP *dr,
                                         const DG_FP *ds, const DG_FP *dt,
                                         const DG_FP *mass, const DG_FP *rx,
                                         const DG_FP *sx, const DG_FP *tx,
                                         const DG_FP *ry, const DG_FP *sy,
                                         const DG_FP *ty, const DG_FP *rz,
                                         const DG_FP *sz, const DG_FP *tz,
                                         const DG_FP *J, const DG_FP *factor,
                                         DG_FP *op1) {
  const DG_FP *dr_mat = &dr[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dt[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *mass_mat = &mass[(*order - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];

  DG_FP D[DG_NP * DG_NP], D_f[DG_NP * DG_NP], D_t[DG_NP * DG_NP];
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      D[ind] = rx[0] * dr_mat[ind] + sx[0] * ds_mat[ind] + tx[0] * dt_mat[ind];
      D_f[ind] = D[ind] * factor[i];
    }
  }

  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, J[0], mass_mat, dg_np, D_f, dg_np, 0.0, D_t, dg_np);
  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, 1.0, D, dg_np, D_t, dg_np, 0.0, op1, dg_np);

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      D[ind] = ry[0] * dr_mat[ind] + sy[0] * ds_mat[ind] + ty[0] * dt_mat[ind];
      D_f[ind] = D[ind] * factor[i];
    }
  }

  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, J[0], mass_mat, dg_np, D_f, dg_np, 0.0, D_t, dg_np);
  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, 1.0, D, dg_np, D_t, dg_np, 1.0, op1, dg_np);

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      D[ind] = rz[0] * dr_mat[ind] + sz[0] * ds_mat[ind] + tz[0] * dt_mat[ind];
      D_f[ind] = D[ind] * factor[i];
    }
  }

  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, J[0], mass_mat, dg_np, D_f, dg_np, 0.0, D_t, dg_np);
  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, 1.0, D, dg_np, D_t, dg_np, 1.0, op1, dg_np);
}
