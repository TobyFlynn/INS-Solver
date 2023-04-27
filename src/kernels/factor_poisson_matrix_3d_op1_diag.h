inline void factor_poisson_matrix_3d_op1_diag(const int *order, const DG_FP *rx,
                                         const DG_FP *sx, const DG_FP *tx,
                                         const DG_FP *ry, const DG_FP *sy,
                                         const DG_FP *ty, const DG_FP *rz,
                                         const DG_FP *sz, const DG_FP *tz,
                                         const DG_FP *J, const DG_FP *factor,
                                         DG_FP *diag) {
  const DG_FP *dr_mat = &dg_Dr_kernel[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &dg_Ds_kernel[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dg_Dt_kernel[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *mass_mat = &dg_Mass_kernel[(*order - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];

  // X

  DG_FP D_t[DG_NP * DG_NP], D[DG_NP * DG_NP];
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      int op_ind = DG_MAT_IND(i, j, dg_np, dg_np);
      D[op_ind] = rx[0] * dr_mat[op_ind] + sx[0] * ds_mat[op_ind] + tx[0] * dt_mat[op_ind];
      D[op_ind] *= factor[i];
    }
  }

  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, 1.0, mass_mat, dg_np, D, dg_np, 0.0, D_t, dg_np);

  for(int i = 0; i < dg_np; i++) {
    diag[i] = 0.0;
    for(int k = 0; k < dg_np; k++) {
      // int a_ind = i * dg_np + k;
      int a_ind = DG_MAT_IND(k, i, dg_np, dg_np);
      // int b_ind = j * dg_np + k;
      int b_ind = DG_MAT_IND(k, i, dg_np, dg_np);
      DG_FP Dx = rx[0] * dr_mat[a_ind] + sx[0] * ds_mat[a_ind] + tx[0] * dt_mat[a_ind];
      diag[i] += Dx * D_t[b_ind];
    }
  }

  // Y
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      int op_ind = DG_MAT_IND(i, j, dg_np, dg_np);
      D[op_ind] = ry[0] * dr_mat[op_ind] + sy[0] * ds_mat[op_ind] + ty[0] * dt_mat[op_ind];
      D[op_ind] *= factor[i];
    }
  }

  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, 1.0, mass_mat, dg_np, D, dg_np, 0.0, D_t, dg_np);

  for(int i = 0; i < dg_np; i++) {
    for(int k = 0; k < dg_np; k++) {
      // int a_ind = i * dg_np + k;
      int a_ind = DG_MAT_IND(k, i, dg_np, dg_np);
      // int b_ind = j * dg_np + k;
      int b_ind = DG_MAT_IND(k, i, dg_np, dg_np);
      DG_FP Dy = ry[0] * dr_mat[a_ind] + sy[0] * ds_mat[a_ind] + ty[0] * dt_mat[a_ind];
      diag[i] += Dy * D_t[b_ind];
    }
  }

  // Z
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      int op_ind = DG_MAT_IND(i, j, dg_np, dg_np);
      D[op_ind] = rz[0] * dr_mat[op_ind] + sz[0] * ds_mat[op_ind] + tz[0] * dt_mat[op_ind];
      D[op_ind] *= factor[i];
    }
  }

  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, 1.0, mass_mat, dg_np, D, dg_np, 0.0, D_t, dg_np);

  for(int i = 0; i < dg_np; i++) {
    for(int k = 0; k < dg_np; k++) {
      // int a_ind = i * dg_np + k;
      int a_ind = DG_MAT_IND(k, i, dg_np, dg_np);
      // int b_ind = j * dg_np + k;
      int b_ind = DG_MAT_IND(k, i, dg_np, dg_np);
      DG_FP Dz = rz[0] * dr_mat[a_ind] + sz[0] * ds_mat[a_ind] + tz[0] * dt_mat[a_ind];
      diag[i] += Dz * D_t[b_ind];
    }
    diag[i] *= J[0];
  }
}
