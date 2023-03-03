inline void factor_poisson_coarse_matrix_3d_op1(const DG_FP *dr,
                    const DG_FP *ds, const DG_FP *dt, const DG_FP *mass,
                    const DG_FP *rx, const DG_FP *sx, const DG_FP *tx,
                    const DG_FP *ry, const DG_FP *sy, const DG_FP *ty,
                    const DG_FP *rz, const DG_FP *sz, const DG_FP *tz,
                    const DG_FP *J, const DG_FP *factor, DG_FP *op1) {
  const DG_FP *dr_mat = dr;
  const DG_FP *ds_mat = ds;
  const DG_FP *dt_mat = dt;
  const DG_FP *mass_mat = mass;

  DG_FP Dx[DG_NP_N1 * DG_NP_N1], Dy[DG_NP_N1 * DG_NP_N1], Dz[DG_NP_N1 * DG_NP_N1];
  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int ind = i + j * DG_NP_N1;
      int ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);
      Dx[ind] = rx[0] * dr_mat[ind] + sx[0] * ds_mat[ind] + tx[0] * dt_mat[ind];
      Dy[ind] = ry[0] * dr_mat[ind] + sy[0] * ds_mat[ind] + ty[0] * dt_mat[ind];
      Dz[ind] = rz[0] * dr_mat[ind] + sz[0] * ds_mat[ind] + tz[0] * dt_mat[ind];
    }
  }

  DG_FP Dx_t[DG_NP_N1 * DG_NP_N1], Dy_t[DG_NP_N1 * DG_NP_N1], Dz_t[DG_NP_N1 * DG_NP_N1];
  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int op_ind = i + j * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);
      Dx_t[op_ind] = 0.0;
      Dy_t[op_ind] = 0.0;
      Dz_t[op_ind] = 0.0;
      for(int k = 0; k < DG_NP_N1; k++) {
        // int a_ind = i + k * DG_NP_N1;
        int a_ind = DG_MAT_IND(i, k, DG_NP_N1, DG_NP_N1);
        // int b_ind = j * DG_NP_N1 + k;
        int b_ind = DG_MAT_IND(k, j, DG_NP_N1, DG_NP_N1);
        Dx_t[op_ind] += mass_mat[a_ind] * Dx[b_ind];
        Dy_t[op_ind] += mass_mat[a_ind] * Dy[b_ind];
        Dz_t[op_ind] += mass_mat[a_ind] * Dz[b_ind];
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
        DG_FP tmp = Dx[a_ind] * Dx_t[b_ind] + Dy[a_ind] * Dy_t[b_ind] + Dz[a_ind] * Dz_t[b_ind];
        op1[op_ind] += factor[k] * tmp;
      }
      op1[op_ind] *= J[0];
    }
  }
}
