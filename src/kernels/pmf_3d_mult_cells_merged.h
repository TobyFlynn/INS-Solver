inline void pmf_3d_mult_cells_merged(const int *p, const DG_FP *emat, const DG_FP *mass,
                            const DG_FP *dr, const DG_FP *ds, const DG_FP *dt,
                            const DG_FP *rx, const DG_FP *sx, const DG_FP *tx,
                            const DG_FP *ry, const DG_FP *sy, const DG_FP *ty,
                            const DG_FP *rz, const DG_FP *sz, const DG_FP *tz,
                            const DG_FP *J, const DG_FP *l_x, const DG_FP *l_y,
                            const DG_FP *l_z, const DG_FP *out_tmp,
                            const DG_FP *in_x, const DG_FP *in_y,
                            const DG_FP *in_z, DG_FP *out) {
  const DG_FP *emat_mat = &emat[(*p - 1) * DG_NUM_FACES * DG_NPF * DG_NP];
  const DG_FP *mass_mat = &mass[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *dr_mat = &dr[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dt[(*p - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 1];

  DG_FP tmp_x[DG_NP], tmp_y[DG_NP], tmp_z[DG_NP];
  op2_in_kernel_gemv(false, dg_np, dg_np, *J, mass_mat, dg_np, in_x, 0.0, tmp_x);
  op2_in_kernel_gemv(false, dg_np, dg_np, *J, mass_mat, dg_np, in_y, 0.0, tmp_y);
  op2_in_kernel_gemv(false, dg_np, dg_np, *J, mass_mat, dg_np, in_z, 0.0, tmp_z);

  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, l_x, 1.0, tmp_x);
  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, l_y, 1.0, tmp_y);
  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, l_z, 1.0, tmp_z);
  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, out_tmp, 0.0, out);

  DG_FP tmp_dr[DG_NP], tmp_ds[DG_NP], tmp_dt[DG_NP];
  for(int n = 0; n < dg_np; n++) {
    tmp_dr[n] = rx[0] * tmp_x[n] + ry[0] * tmp_y[n] + rz[0] * tmp_z[n];
    tmp_ds[n] = sx[0] * tmp_x[n] + sy[0] * tmp_y[n] + sz[0] * tmp_z[n];
    tmp_dt[n] = tx[0] * tmp_x[n] + ty[0] * tmp_y[n] + tz[0] * tmp_z[n];
  }

  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, dr_mat, dg_np, tmp_dr, 1.0, out);
  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, ds_mat, dg_np, tmp_ds, 1.0, out);
  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, dt_mat, dg_np, tmp_dt, 1.0, out);
}
