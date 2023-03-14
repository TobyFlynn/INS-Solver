inline void pmf_2d_mult_cells(const int *p, const DG_FP *dr, const DG_FP *ds,
                      const DG_FP *rx, const DG_FP *sx, const DG_FP *ry,
                      const DG_FP *sy, const DG_FP *l_x, const DG_FP *l_y,
                      const DG_FP *in_x, const DG_FP *in_y, DG_FP *out) {
  const DG_FP *dr_mat = &dr[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*p - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  DG_FP tmp_dr[DG_NP], tmp_ds[DG_NP];
  for(int n = 0; n < dg_np; n++) {
    tmp_dr[n] = rx[n] * (in_x[n] + l_x[n]) + ry[n] * (in_y[n] + l_y[n]);
    tmp_ds[n] = sx[n] * (in_x[n] + l_x[n]) + sy[n] * (in_y[n] + l_y[n]);
  }

  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, dr_mat, dg_np, tmp_dr, 1.0, out);
  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, ds_mat, dg_np, tmp_ds, 1.0, out);
}
