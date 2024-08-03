inline void fvmf_2d_op1_diag(const DG_FP *geof, const DG_FP *factor, DG_FP *diag) {
  const DG_FP *dr_mat = &dg_Dr_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &dg_Ds_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mass_mat = &dg_Mass_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];

  // X

  DG_FP D_t[DG_NP * DG_NP], D[DG_NP * DG_NP];
  const DG_FP rx = geof[RX_IND];
  const DG_FP sx = geof[SX_IND];
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int op_ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      D[op_ind] = rx * dr_mat[op_ind] + sx * ds_mat[op_ind];
      D[op_ind] *= factor[i];
    }
  }

  op2_in_kernel_gemm(false, false, DG_NP, DG_NP, DG_NP, 1.0, mass_mat, DG_NP, D, DG_NP, 0.0, D_t, DG_NP);

  for(int i = 0; i < DG_NP; i++) {
    diag[i] = 0.0;
    for(int k = 0; k < DG_NP; k++) {
      // int a_ind = i * DG_NP + k;
      int a_ind = DG_MAT_IND(k, i, DG_NP, DG_NP);
      // int b_ind = j * DG_NP + k;
      int b_ind = DG_MAT_IND(k, i, DG_NP, DG_NP);
      DG_FP Dx = rx * dr_mat[a_ind] + sx * ds_mat[a_ind];
      diag[i] += Dx * D_t[b_ind];
    }
  }

  // Y
  const DG_FP ry = geof[RY_IND];
  const DG_FP sy = geof[SY_IND];
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int op_ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      D[op_ind] = ry * dr_mat[op_ind] + sy * ds_mat[op_ind];
      D[op_ind] *= factor[i];
    }
  }

  op2_in_kernel_gemm(false, false, DG_NP, DG_NP, DG_NP, 1.0, mass_mat, DG_NP, D, DG_NP, 0.0, D_t, DG_NP);

  const DG_FP J = geof[J_IND];
  for(int i = 0; i < DG_NP; i++) {
    for(int k = 0; k < DG_NP; k++) {
      // int a_ind = i * DG_NP + k;
      int a_ind = DG_MAT_IND(k, i, DG_NP, DG_NP);
      // int b_ind = j * DG_NP + k;
      int b_ind = DG_MAT_IND(k, i, DG_NP, DG_NP);
      DG_FP Dy = ry * dr_mat[a_ind] + sy * ds_mat[a_ind];
      diag[i] += Dy * D_t[b_ind];
    }
    diag[i] *= J;
  }
}