inline void pmf_3d_mult_complete_flux(const int *p, const DG_FP *emat, const DG_FP *dr,
                             const DG_FP *ds, const DG_FP *dt, const DG_FP *mass,
                             const int *faceNums, const int *fmaskF,
                             const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                             const DG_FP *fscale, const DG_FP *sJ,
                             const DG_FP *rx, const DG_FP *sx, const DG_FP *tx,
                             const DG_FP *ry, const DG_FP *sy, const DG_FP *ty,
                             const DG_FP *rz, const DG_FP *sz, const DG_FP *tz,
                             const DG_FP *J, const DG_FP **in, const DG_FP **in_x,
                             const DG_FP **in_y, const DG_FP **in_z, DG_FP *out) {
  const int dg_np  = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 1];
  const int *fmask  = &FMASK[(*p - 1) * 4 * DG_NPF];
  const DG_FP *emat_mat = &emat[(*p - 1) * DG_NUM_FACES * DG_NPF * DG_NP];
  const DG_FP *dr_mat = &dr[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dt[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *mass_mat = &mass[(*p - 1) * DG_NP * DG_NP];

  DG_FP tmp_0[DG_NPF * DG_NUM_FACES], tmp_1[DG_NPF * DG_NUM_FACES], tmp_2[DG_NPF * DG_NUM_FACES], tmp_3[DG_NPF * DG_NUM_FACES];
  for(int i = 0; i < DG_NUM_FACES; i++) {
    const DG_FP gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(fscale[i * 2], fscale[i * 2 + 1]);

    const int findL = faceNums[2 * i] * dg_npf;
    const int *fmaskL = &fmask[faceNums[2 * i] * dg_npf];
    const int *fmaskR = &fmaskF[i * dg_npf];

    for(int j = 0; j < dg_npf; j++) {
      const DG_FP diffL_u = in[0][fmaskL[j]] - in[i + 1][fmaskR[j]];
      const DG_FP diffL_u_x = nx[i] * (in_x[i + 1][fmaskR[j]] + in_x[0][fmaskL[j]]);
      const DG_FP diffL_u_y = ny[i] * (in_y[i + 1][fmaskR[j]] + in_y[0][fmaskL[j]]);
      const DG_FP diffL_u_z = nz[i] * (in_z[i + 1][fmaskR[j]] + in_z[0][fmaskL[j]]);
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
  op2_in_kernel_gemv(false, dg_np, dg_np, J[0], mass_mat, dg_np, in_x[0], 1.0, tmp_4);
  op2_in_kernel_gemv(false, dg_np, dg_np, J[0], mass_mat, dg_np, in_y[0], 1.0, tmp_5);
  op2_in_kernel_gemv(false, dg_np, dg_np, J[0], mass_mat, dg_np, in_z[0], 1.0, tmp_6);

  DG_FP tmp_dr[DG_NP], tmp_ds[DG_NP], tmp_dt[DG_NP];
  for(int n = 0; n < dg_np; n++) {
    tmp_dr[n] = rx[0] * tmp_4[n] + ry[0] * tmp_5[n] + rz[0] * tmp_6[n];
    tmp_ds[n] = sx[0] * tmp_4[n] + sy[0] * tmp_5[n] + sz[0] * tmp_6[n];
    tmp_dt[n] = tx[0] * tmp_4[n] + ty[0] * tmp_5[n] + tz[0] * tmp_6[n];
  }

  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, dr_mat, dg_np, tmp_dr, 1.0, out);
  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, ds_mat, dg_np, tmp_ds, 1.0, out);
  op2_in_kernel_gemv(true, dg_np, dg_np, 1.0, dt_mat, dg_np, tmp_dt, 1.0, out);
}
