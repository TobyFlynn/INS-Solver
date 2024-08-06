inline void fvmf_3d_op1(const DG_FP *geof, const DG_FP *factor, DG_FP *op1) {
  const DG_FP *dr_mat = &dg_Dr_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &dg_Ds_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dg_Dt_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mass_mat = &dg_Mass_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];

  DG_FP D[DG_NP * DG_NP], D_f[DG_NP * DG_NP], D_t[DG_NP * DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      // int ind = i + j * DG_NP;
      int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      D[ind] = geof[RX_IND] * dr_mat[ind] + geof[SX_IND] * ds_mat[ind] + geof[TX_IND] * dt_mat[ind];
      D_f[ind] = D[ind] * factor[i];
    }
  }

  op2_in_kernel_gemm(false, false, DG_NP, DG_NP, DG_NP, geof[J_IND], mass_mat, DG_NP, D_f, DG_NP, 0.0, D_t, DG_NP);
  op2_in_kernel_gemm(false, false, DG_NP, DG_NP, DG_NP, 1.0, D, DG_NP, D_t, DG_NP, 0.0, op1, DG_NP);

  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      // int ind = i + j * DG_NP;
      int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      D[ind] = geof[RY_IND] * dr_mat[ind] + geof[SY_IND] * ds_mat[ind] + geof[TY_IND] * dt_mat[ind];
      D_f[ind] = D[ind] * factor[i];
    }
  }

  op2_in_kernel_gemm(false, false, DG_NP, DG_NP, DG_NP, geof[J_IND], mass_mat, DG_NP, D_f, DG_NP, 0.0, D_t, DG_NP);
  op2_in_kernel_gemm(false, false, DG_NP, DG_NP, DG_NP, 1.0, D, DG_NP, D_t, DG_NP, 1.0, op1, DG_NP);

  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      // int ind = i + j * DG_NP;
      int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      D[ind] = geof[RZ_IND] * dr_mat[ind] + geof[SZ_IND] * ds_mat[ind] + geof[TZ_IND] * dt_mat[ind];
      D_f[ind] = D[ind] * factor[i];
    }
  }

  op2_in_kernel_gemm(false, false, DG_NP, DG_NP, DG_NP, geof[J_IND], mass_mat, DG_NP, D_f, DG_NP, 0.0, D_t, DG_NP);
  op2_in_kernel_gemm(false, false, DG_NP, DG_NP, DG_NP, 1.0, D, DG_NP, D_t, DG_NP, 1.0, op1, DG_NP);
}