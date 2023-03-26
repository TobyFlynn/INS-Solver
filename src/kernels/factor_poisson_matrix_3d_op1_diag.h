inline void factor_poisson_matrix_3d_op1_diag(const int *order, const DG_FP *dr,
                                         const DG_FP *ds, const DG_FP *dt,
                                         const DG_FP *mass, const DG_FP *rx,
                                         const DG_FP *sx, const DG_FP *tx,
                                         const DG_FP *ry, const DG_FP *sy,
                                         const DG_FP *ty, const DG_FP *rz,
                                         const DG_FP *sz, const DG_FP *tz,
                                         const DG_FP *J, const DG_FP *factor,
                                         DG_FP *diag) {
  const DG_FP *dr_mat = &dr[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dt[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *mass_mat = &mass[(*order - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];

  DG_FP D_t[DG_NP * DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      // int op_ind = i + j * DG_NP;
      int op_ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      D_t[op_ind] = 0.0;
      for(int k = 0; k < DG_NP; k++) {
        // int a_ind = i + k * DG_NP;
        int a_ind = DG_MAT_IND(i, k, DG_NP, DG_NP);
        // int b_ind = j * DG_NP + k;
        int b_ind = DG_MAT_IND(k, j, DG_NP, DG_NP);
        DG_FP Dx = rx[0] * dr_mat[b_ind] + sx[0] * ds_mat[b_ind] + tx[0] * dt_mat[b_ind];
        D_t[op_ind] += mass_mat[a_ind] * factor[k] * Dx;
      }
    }
  }

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

  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      // int op_ind = i + j * DG_NP;
      int op_ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      D_t[op_ind] = 0.0;
      for(int k = 0; k < DG_NP; k++) {
        // int a_ind = i + k * DG_NP;
        int a_ind = DG_MAT_IND(i, k, DG_NP, DG_NP);
        // int b_ind = j * DG_NP + k;
        int b_ind = DG_MAT_IND(k, j, DG_NP, DG_NP);
        DG_FP Dy = ry[0] * dr_mat[b_ind] + sy[0] * ds_mat[b_ind] + ty[0] * dt_mat[b_ind];
        D_t[op_ind] += mass_mat[a_ind] * factor[k] * Dy;
      }
    }
  }

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

  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      // int op_ind = i + j * DG_NP;
      int op_ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      D_t[op_ind] = 0.0;
      for(int k = 0; k < DG_NP; k++) {
        // int a_ind = i + k * DG_NP;
        int a_ind = DG_MAT_IND(i, k, DG_NP, DG_NP);
        // int b_ind = j * DG_NP + k;
        int b_ind = DG_MAT_IND(k, j, DG_NP, DG_NP);
        DG_FP Dz = rz[0] * dr_mat[b_ind] + sz[0] * ds_mat[b_ind] + tz[0] * dt_mat[b_ind];
        D_t[op_ind] += mass_mat[a_ind] * factor[k] * Dz;
      }
    }
  }

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
