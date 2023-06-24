inline void pmf_3d_mult_complete_flux(const int *p, const int *faceNums,
                            const int *fmaskF, const DG_FP *nx, const DG_FP *ny,
                            const DG_FP *nz, const DG_FP *fscale,
                            const DG_FP *sJ, const DG_FP *geof, const DG_FP *in,
                            const DG_FP **in_p, const DG_FP *in_x,
                            const DG_FP *in_y, const DG_FP *in_z,
                            const DG_FP **in_x_p, const DG_FP **in_y_p,
                            const DG_FP **in_z_p, DG_FP *out) {
  const int dg_np  = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 1];
  const int *fmask  = &FMASK[(*p - 1) * 4 * DG_NPF];
  const DG_FP *emat_mat = &dg_Emat_kernel[(*p - 1) * DG_NUM_FACES * DG_NPF * DG_NP];
  const DG_FP *dr_mat = &dg_Dr_kernel[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &dg_Ds_kernel[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dg_Dt_kernel[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *mass_mat = &dg_Mass_kernel[(*p - 1) * DG_NP * DG_NP];

  DG_FP tmp_0[DG_NPF * DG_NUM_FACES], tmp_1[DG_NPF * DG_NUM_FACES], tmp_2[DG_NPF * DG_NUM_FACES], tmp_3[DG_NPF * DG_NUM_FACES];
  for(int i = 0; i < DG_NUM_FACES; i++) {
    const DG_FP gtau = 2.0 * (*p + 1) * (*p + 2) * fmax(fscale[i * 2], fscale[i * 2 + 1]);

    const int findL = faceNums[2 * i] * dg_npf;
    const int *fmaskL = &fmask[faceNums[2 * i] * dg_npf];

    for(int j = 0; j < dg_npf; j++) {
      const int fmaskIndL = fmaskL[j];
      const int fmaskIndR = fmaskF[i * dg_npf + j];
      const DG_FP diffL_u = in[fmaskIndL] - in_p[i][fmaskIndR];
      const DG_FP diffL_u_x = nx[i] * (in_x_p[i][fmaskIndR] + in_x[fmaskIndL]);
      const DG_FP diffL_u_y = ny[i] * (in_y_p[i][fmaskIndR] + in_y[fmaskIndL]);
      const DG_FP diffL_u_z = nz[i] * (in_z_p[i][fmaskIndR] + in_z[fmaskIndL]);
      const DG_FP diffL_u_grad = diffL_u_x + diffL_u_y + diffL_u_z;

      const int indL = findL + j;
      tmp_0[indL] = 0.5 * sJ[i] * (gtau * diffL_u - diffL_u_grad);
      const DG_FP l_tmpL = 0.5 * sJ[i] * -diffL_u;
      tmp_1[indL] = nx[i] * l_tmpL;
      tmp_2[indL] = ny[i] * l_tmpL;
      tmp_3[indL] = nz[i] * l_tmpL;
    }
  }

  DG_FP tmp_4[DG_NP], tmp_5[DG_NP], tmp_6[DG_NP];
  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, tmp_0, 0.0, out);
  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, tmp_1, 0.0, tmp_4);
  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, tmp_2, 0.0, tmp_5);
  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, tmp_3, 0.0, tmp_6);

  // MASS
  op2_in_kernel_gemv(false, dg_np, dg_np, geof[J_IND], mass_mat, dg_np, in_x, 1.0, tmp_4);
  op2_in_kernel_gemv(false, dg_np, dg_np, geof[J_IND], mass_mat, dg_np, in_y, 1.0, tmp_5);
  op2_in_kernel_gemv(false, dg_np, dg_np, geof[J_IND], mass_mat, dg_np, in_z, 1.0, tmp_6);

  DG_FP tmp_dr[DG_NP], tmp_ds[DG_NP], tmp_dt[DG_NP];
  for(int n = 0; n < dg_np; n++) {
    tmp_dr[n] = geof[RX_IND] * tmp_4[n] + geof[RY_IND] * tmp_5[n] + geof[RZ_IND] * tmp_6[n];
    tmp_ds[n] = geof[SX_IND] * tmp_4[n] + geof[SY_IND] * tmp_5[n] + geof[SZ_IND] * tmp_6[n];
    tmp_dt[n] = geof[TX_IND] * tmp_4[n] + geof[TY_IND] * tmp_5[n] + geof[TZ_IND] * tmp_6[n];
  }

  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, dr_mat, dg_np, tmp_dr, 1.0, out);
  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, ds_mat, dg_np, tmp_ds, 1.0, out);
  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, dt_mat, dg_np, tmp_dt, 1.0, out);
}
