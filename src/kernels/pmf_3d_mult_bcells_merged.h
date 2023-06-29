inline void pmf_3d_mult_bcells_merged(const int *p, const int *bcell,
                          const DG_FP *geof, const DG_FP *l_x, const DG_FP *l_y,
                          const DG_FP *l_z, const DG_FP *out_tmp,
                          const DG_FP *in_x, const DG_FP *in_y,
                          const DG_FP *in_z, DG_FP *out) {
  if(*bcell == 0) return;
  const DG_FP *emat_mat = &dg_Emat_kernel[(*p - 1) * DG_NUM_FACES * DG_NPF * DG_NP];
  const DG_FP *mass_mat = &dg_Mass_kernel[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *dr_mat = &dg_Dr_kernel[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &dg_Ds_kernel[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dg_Dt_kernel[(*p - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 1];

  DG_FP tmp_x[DG_NP], tmp_y[DG_NP], tmp_z[DG_NP];
  op2_in_kernel_gemv(false, dg_np, dg_np, geof[J_IND], mass_mat, dg_np, in_x, 0.0, tmp_x);
  op2_in_kernel_gemv(false, dg_np, dg_np, geof[J_IND], mass_mat, dg_np, in_y, 0.0, tmp_y);
  op2_in_kernel_gemv(false, dg_np, dg_np, geof[J_IND], mass_mat, dg_np, in_z, 0.0, tmp_z);

  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, l_x, 1.0, tmp_x);
  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, l_y, 1.0, tmp_y);
  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, l_z, 1.0, tmp_z);
  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, out_tmp, 0.0, out);

  DG_FP tmp_dr[DG_NP], tmp_ds[DG_NP], tmp_dt[DG_NP];
  for(int n = 0; n < dg_np; n++) {
    tmp_dr[n] = geof[RX_IND] * tmp_x[n] + geof[RY_IND] * tmp_y[n] + geof[RZ_IND] * tmp_z[n];
    tmp_ds[n] = geof[SX_IND] * tmp_x[n] + geof[SY_IND] * tmp_y[n] + geof[SZ_IND] * tmp_z[n];
    tmp_dt[n] = geof[TX_IND] * tmp_x[n] + geof[TY_IND] * tmp_y[n] + geof[TZ_IND] * tmp_z[n];
  }

  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, dr_mat, dg_np, tmp_dr, 1.0, out);
  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, ds_mat, dg_np, tmp_ds, 1.0, out);
  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, dt_mat, dg_np, tmp_dt, 1.0, out);
}
