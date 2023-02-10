inline void poisson_matrix_3d_op1(const int *order, const double *dr,
                                  const double *ds, const double *dt,
                                  const double *mass, const double *rx,
                                  const double *sx, const double *tx,
                                  const double *ry, const double *sy,
                                  const double *ty, const double *rz,
                                  const double *sz, const double *tz,
                                  const double *J, double *op1) {
  const double *dr_mat = &dr[(*order - 1) * DG_NP * DG_NP];
  const double *ds_mat = &ds[(*order - 1) * DG_NP * DG_NP];
  const double *dt_mat = &dt[(*order - 1) * DG_NP * DG_NP];
  const double *mass_mat = &mass[(*order - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];

  double Dx[DG_NP * DG_NP], Dy[DG_NP * DG_NP], Dz[DG_NP * DG_NP];
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      int ind = i + j * dg_np;
      Dx[ind] = rx[0] * dr_mat[ind] + sx[0] * ds_mat[ind] + tx[0] * dt_mat[ind];
      Dy[ind] = ry[0] * dr_mat[ind] + sy[0] * ds_mat[ind] + ty[0] * dt_mat[ind];
      Dz[ind] = rz[0] * dr_mat[ind] + sz[0] * ds_mat[ind] + tz[0] * dt_mat[ind];
    }
  }

  double Dx_t[DG_NP * DG_NP], Dy_t[DG_NP * DG_NP], Dz_t[DG_NP * DG_NP];
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      int op_ind = i + j * dg_np;
      Dx_t[op_ind] = 0.0;
      Dy_t[op_ind] = 0.0;
      Dz_t[op_ind] = 0.0;
      for(int k = 0; k < dg_np; k++) {
        int a_ind = i + k * dg_np;
        int b_ind = j * dg_np + k;
        Dx_t[op_ind] += mass_mat[a_ind] * Dx[b_ind];
        Dy_t[op_ind] += mass_mat[a_ind] * Dy[b_ind];
        Dz_t[op_ind] += mass_mat[a_ind] * Dz[b_ind];
      }
    }
  }

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      int op_ind = i + j * dg_np;
      op1[op_ind] = 0.0;
      for(int k = 0; k < dg_np; k++) {
        int a_ind = i * dg_np + k;
        int b_ind = j * dg_np + k;
        op1[op_ind] += Dx[a_ind] * Dx_t[b_ind] + Dy[a_ind] * Dy_t[b_ind] + Dz[a_ind] * Dz_t[b_ind];
      }
      op1[op_ind] *= J[0];
    }
  }
}
