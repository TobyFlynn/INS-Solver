inline void fvmf_2d_op1(const DG_FP *geof, const DG_FP *factor, DG_FP *op1) {
  const DG_FP *dr_mat = &dg_Dr_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &dg_Ds_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mass_mat = &dg_Mass_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];

  DG_FP Dx[DG_NP * DG_NP], Dy[DG_NP * DG_NP];
  const DG_FP rx = geof[RX_IND];
  const DG_FP sx = geof[SX_IND];
  const DG_FP ry = geof[RY_IND];
  const DG_FP sy = geof[SY_IND];
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      // int ind = i + j * DG_NP;
      int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      Dx[ind] = rx * dr_mat[ind] + sx * ds_mat[ind];
      Dy[ind] = ry * dr_mat[ind] + sy * ds_mat[ind];
    }
  }

  DG_FP Dx_t[DG_NP * DG_NP], Dy_t[DG_NP * DG_NP];
  op2_in_kernel_gemm(false, false, DG_NP, DG_NP, DG_NP, 1.0, mass_mat, DG_NP, Dx, DG_NP, 0.0, Dx_t, DG_NP);
  op2_in_kernel_gemm(false, false, DG_NP, DG_NP, DG_NP, 1.0, mass_mat, DG_NP, Dy, DG_NP, 0.0, Dy_t, DG_NP);

  const DG_FP J = geof[J_IND];
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      // int op_ind = i + j * DG_NP;
      int op_ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      op1[op_ind] = 0.0;
      for(int k = 0; k < DG_NP; k++) {
        // int a_ind = i * DG_NP + k;
        int a_ind = DG_MAT_IND(k, i, DG_NP, DG_NP);
        // int b_ind = j * DG_NP + k;
        int b_ind = DG_MAT_IND(k, j, DG_NP, DG_NP);
        DG_FP tmp = Dx[a_ind] * Dx_t[b_ind] + Dy[a_ind] * Dy_t[b_ind];
        op1[op_ind] += factor[k] * tmp;
      }
      op1[op_ind] *= J;
    }
  }
}