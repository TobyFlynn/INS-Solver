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

  DG_FP Dx[DG_NP * DG_NP], Dy[DG_NP * DG_NP], Dz[DG_NP * DG_NP];
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      Dx[ind] = rx[0] * dr_mat[ind] + sx[0] * ds_mat[ind] + tx[0] * dt_mat[ind];
      Dy[ind] = ry[0] * dr_mat[ind] + sy[0] * ds_mat[ind] + ty[0] * dt_mat[ind];
      Dz[ind] = rz[0] * dr_mat[ind] + sz[0] * ds_mat[ind] + tz[0] * dt_mat[ind];
    }
  }

  DG_FP Dx_t[DG_NP * DG_NP], Dy_t[DG_NP * DG_NP], Dz_t[DG_NP * DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      // int op_ind = i + j * DG_NP;
      int op_ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      Dx_t[op_ind] = 0.0;
      Dy_t[op_ind] = 0.0;
      Dz_t[op_ind] = 0.0;
      for(int k = 0; k < DG_NP; k++) {
        // int a_ind = i + k * DG_NP;
        int a_ind = DG_MAT_IND(i, k, DG_NP, DG_NP);
        // int b_ind = j * DG_NP + k;
        int b_ind = DG_MAT_IND(k, j, DG_NP, DG_NP);
        Dx_t[op_ind] += mass_mat[a_ind] * factor[k] * Dx[b_ind];
        Dy_t[op_ind] += mass_mat[a_ind] * factor[k] * Dy[b_ind];
        Dz_t[op_ind] += mass_mat[a_ind] * factor[k] * Dz[b_ind];
      }
    }
  }

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int op_ind = i + j * dg_np;
      int op_ind = DG_MAT_IND(i, j, dg_np, dg_np);
      op1[op_ind] = 0.0;
      for(int k = 0; k < dg_np; k++) {
        // int a_ind = i * dg_np + k;
        int a_ind = DG_MAT_IND(k, i, dg_np, dg_np);
        // int b_ind = j * dg_np + k;
        int b_ind = DG_MAT_IND(k, j, dg_np, dg_np);
        DG_FP tmp = Dx[a_ind] * Dx_t[b_ind] + Dy[a_ind] * Dy_t[b_ind] + Dz[a_ind] * Dz_t[b_ind];
        op1[op_ind] += factor[k] * tmp;
      }
      op1[op_ind] *= J[0];
    }
  }
/*
  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, 1.0, mass_mat, dg_np, Dx, dg_np, 0.0, Dx_t, dg_np);
  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, 1.0, mass_mat, dg_np, Dy, dg_np, 0.0, Dy_t, dg_np);
  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, 1.0, mass_mat, dg_np, Dz, dg_np, 0.0, Dz_t, dg_np);

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int op_ind = i + j * dg_np;
      int op_ind = DG_MAT_IND(i, j, dg_np, dg_np);
      op1[op_ind] = 0.0;
      for(int k = 0; k < dg_np; k++) {
        // int a_ind = i * dg_np + k;
        int a_ind = DG_MAT_IND(k, i, dg_np, dg_np);
        // int b_ind = j * dg_np + k;
        int b_ind = DG_MAT_IND(k, j, dg_np, dg_np);
        DG_FP tmp = Dx[a_ind] * Dx_t[b_ind] + Dy[a_ind] * Dy_t[b_ind] + Dz[a_ind] * Dz_t[b_ind];
        op1[op_ind] += factor[k] * tmp;
      }
      op1[op_ind] *= J[0];
    }
  }
*/
}
