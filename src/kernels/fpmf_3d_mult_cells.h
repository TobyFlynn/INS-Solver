inline void fpmf_3d_mult_cells(const int *p, const DG_FP *dr, const DG_FP *ds,
                               const DG_FP *dt, const DG_FP *rx, const DG_FP *sx,
                               const DG_FP *tx, const DG_FP *ry, const DG_FP *sy,
                               const DG_FP *ty, const DG_FP *rz, const DG_FP *sz,
                               const DG_FP *tz, const DG_FP *factor,
                               const DG_FP *l_x, const DG_FP *l_y,
                               const DG_FP *l_z, const DG_FP *in_x,
                               const DG_FP *in_y, const DG_FP *in_z, DG_FP *out) {
  const DG_FP *dr_mat = &dr[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dt[(*p - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  DG_FP tmp_dr[DG_NP], tmp_ds[DG_NP], tmp_dt[DG_NP];
  for(int n = 0; n < dg_np; n++) {
    tmp_dr[n] = rx[0] * (factor[n] * in_x[n] + l_x[n]) + ry[0] * (factor[n] * in_y[n] + l_y[n]) + rz[0] * (factor[n] * in_z[n] + l_z[n]);
    tmp_ds[n] = sx[0] * (factor[n] * in_x[n] + l_x[n]) + sy[0] * (factor[n] * in_y[n] + l_y[n]) + sz[0] * (factor[n] * in_z[n] + l_z[n]);
    tmp_dt[n] = tx[0] * (factor[n] * in_x[n] + l_x[n]) + ty[0] * (factor[n] * in_y[n] + l_y[n]) + tz[0] * (factor[n] * in_z[n] + l_z[n]);
  }

  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, dr_mat, dg_np, tmp_dr, 1.0, out);
  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, ds_mat, dg_np, tmp_ds, 1.0, out);
  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, dt_mat, dg_np, tmp_dt, 1.0, out);

/*
  for(int n = 0; n < dg_np; n++) {
    for(int m = 0; m < dg_np; m++) {
      // int ind = m * dg_np + n;
      int ind = DG_MAT_IND(n, m, dg_np, dg_np);
      out[m] += dr_mat[ind] * (rx[0] * (factor[n] * in_x[n] + l_x[n]) + ry[0] * (factor[n] * in_y[n] + l_y[n]) + rz[0] * (factor[n] * in_z[n] + l_z[n]));
      out[m] += ds_mat[ind] * (sx[0] * (factor[n] * in_x[n] + l_x[n]) + sy[0] * (factor[n] * in_y[n] + l_y[n]) + sz[0] * (factor[n] * in_z[n] + l_z[n]));
      out[m] += dt_mat[ind] * (tx[0] * (factor[n] * in_x[n] + l_x[n]) + ty[0] * (factor[n] * in_y[n] + l_y[n]) + tz[0] * (factor[n] * in_z[n] + l_z[n]));
    }
  }
*/
}
