inline void poisson_matrix_free_3d_mult_cells_new(const int *p, const DG_FP *dr,
                            const DG_FP *ds, const DG_FP *dt,
                            const DG_FP *rx, const DG_FP *sx, const DG_FP *tx,
                            const DG_FP *ry, const DG_FP *sy, const DG_FP *ty,
                            const DG_FP *rz, const DG_FP *sz, const DG_FP *tz,
                            const DG_FP *J, const DG_FP *in_x, const DG_FP *in_y,
                            const DG_FP *in_z, const DG_FP *l_x, const DG_FP *l_y,
                            const DG_FP *l_z, DG_FP *out) {
  const DG_FP *dr_mat = &dr[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dt[(*p - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  for(int m = 0; m < dg_np; m++) {
    out[m] = 0.0;
    for(int n = 0; n < dg_np; n++) {
      // int ind = m * dg_np + n;
      int ind = DG_MAT_IND(m, n, dg_np, dg_np);
      out[m] += dr_mat[ind] * (rx[0] * (in_x[n] + l_x[n]) + ry[0] * (in_y[n] + l_y[n]) + rz[0] * (in_z[n] + l_z[n]));
      out[m] += ds_mat[ind] * (sx[0] * (in_x[n] + l_x[n]) + sy[0] * (in_y[n] + l_y[n]) + sz[0] * (in_z[n] + l_z[n]));
      out[m] += dt_mat[ind] * (tx[0] * (in_x[n] + l_x[n]) + ty[0] * (in_y[n] + l_y[n]) + tz[0] * (in_z[n] + l_z[n]));
    }
  }
}
