inline void factor_poisson_coarse_matrix_2d_op2(const int *faceNum,
                    const bool *reverse, const DG_FP *nx, const DG_FP *ny,
                    const DG_FP *fscale, const DG_FP *sJ, const DG_FP **geof,
                    const DG_FP **factor, DG_FP *op1L, DG_FP *op1R, DG_FP *op2L,
                    DG_FP *op2R) {
  const DG_FP *dr_mat = dg_Dr_kernel;
  const DG_FP *ds_mat = dg_Ds_kernel;
  const DG_FP *mmF0_mat = dg_MM_F0_kernel;
  const DG_FP *mmF1_mat = dg_MM_F1_kernel;
  const DG_FP *mmF2_mat = dg_MM_F2_kernel;

  const DG_FP tau_order = 1.0; // (DG_FP) DG_ORDER;

  const DG_FP *mmFL, *mmFR;
  if(faceNum[0] == 0)
    mmFL = mmF0_mat;
  else if(faceNum[0] == 1)
    mmFL = mmF1_mat;
  else
    mmFL = mmF2_mat;

  if(faceNum[1] == 0)
    mmFR = mmF0_mat;
  else if(faceNum[1] == 1)
    mmFR = mmF1_mat;
  else
    mmFR = mmF2_mat;

  const int faceNumL = faceNum[0];
  const int faceNumR = faceNum[1];
  const int findL = faceNumL * DG_NPF_N1;
  const int findR = faceNumR * DG_NPF_N1;
  const int *fmask  = FMASK;
  const int *fmaskL = &fmask[faceNumL * DG_NPF_N1];
  const int *fmaskR = &fmask[faceNumR * DG_NPF_N1];

  DG_FP DL[DG_NP_N1 * DG_NP_N1], DR[DG_NP_N1 * DG_NP_N1];
  const DG_FP r_fact_0 = nx[0] * geof[0][RX_IND] + ny[0] * geof[0][RY_IND];
  const DG_FP s_fact_0 = nx[0] * geof[0][SX_IND] + ny[0] * geof[0][SY_IND];
  const DG_FP r_fact_1 = nx[1] * geof[1][RX_IND] + ny[1] * geof[1][RY_IND];
  const DG_FP s_fact_1 = nx[1] * geof[1][SX_IND] + ny[1] * geof[1][SY_IND];
  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int ind = i + j * DG_NP_N1;
      int ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);

      DL[ind] = r_fact_0 * dr_mat[ind] + s_fact_0 * ds_mat[ind];
      DL[ind] *= factor[0][i];

      DR[ind] = r_fact_1 * dr_mat[ind] + s_fact_1 * ds_mat[ind];
      DR[ind] *= factor[1][i];
    }
  }

  DG_FP gtau = 0.0;
  for(int i = 0; i < DG_NPF_N1; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = *reverse ? fmaskR[DG_NPF_N1 - i - 1] : fmaskR[i];
    DG_FP tmp = 2.0 * (tau_order + 1) * (tau_order + 2) * fmax(fscale[0] * factor[0][fmaskL_ind], fscale[1] * factor[1][fmaskR_ind]);
    gtau = fmax(gtau, tmp);
  }

  op2_in_kernel_gemm(false, false, DG_NP_N1, DG_NP_N1, DG_NP_N1, -0.5 * sJ[0], mmFL, DG_NP_N1, DL, DG_NP_N1, 1.0, op1L, DG_NP_N1);
  op2_in_kernel_gemm(true, false, DG_NP_N1, DG_NP_N1, DG_NP_N1, -0.5 * sJ[0], DL, DG_NP_N1, mmFL, DG_NP_N1, 1.0, op1L, DG_NP_N1);

  const DG_FP tmp_constL = 0.5 * sJ[0] * gtau;
  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int op_ind = i + j * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);
      op1L[op_ind] += tmp_constL * mmFL[op_ind];
    }
  }

  for(int i = 0; i < DG_NP_N1 * DG_NP_N1; i++) {
    op2L[i] = 0.0;
  }

  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NPF_N1; j++) {
      const int fmaskL_ind = fmaskL[j];
      const int fmaskR_ind = *reverse ? fmaskR[DG_NPF_N1 - j - 1] : fmaskR[j];
      // int op_ind = i + fmaskR_corrected[j] * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, fmaskR_ind, DG_NP_N1, DG_NP_N1);
      // int find = i + fmaskL[j] * DG_NP_N1;
      int find = DG_MAT_IND(i, fmaskL_ind, DG_NP_N1, DG_NP_N1);
      op2L[op_ind] -= 0.5 * gtau * sJ[0] * mmFL[find];
    }
  }

  for(int i = 0; i < DG_NPF_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      const int fmaskL_ind_i = fmaskL[i];
      // int op_ind = fmaskL[i] + j * DG_NP_N1;
      int op_ind = DG_MAT_IND(fmaskL_ind_i, j, DG_NP_N1, DG_NP_N1);
      for(int k = 0; k < DG_NPF_N1; k++) {
        const int fmaskL_ind_k = fmaskL[k];
        const int fmaskR_ind_k = *reverse ? fmaskR[DG_NPF_N1 - k - 1] : fmaskR[k];
        // int a_ind = fmaskL[i] + fmaskL[k] * DG_NP_N1;
        int a_ind = DG_MAT_IND(fmaskL_ind_i, fmaskL_ind_k, DG_NP_N1, DG_NP_N1);
        // int b_ind = j * DG_NP_N1 + fmaskR_corrected[k];
        int b_ind = DG_MAT_IND(fmaskR_ind_k, j, DG_NP_N1, DG_NP_N1);
        op2L[op_ind] -= 0.5 * sJ[0] * mmFL[a_ind] * -DR[b_ind];
      }
    }
  }

  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NPF_N1; j++) {
      const int fmaskR_ind = *reverse ? fmaskR[DG_NPF_N1 - j - 1] : fmaskR[j];
      // int op_ind = i + fmaskR_corrected[j] * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, fmaskR_ind, DG_NP_N1, DG_NP_N1);
      for(int k = 0; k < DG_NP_N1; k++) {
        const int fmaskL_ind = fmaskL[j];
        // int a_ind = i * DG_NP_N1 + k;
        int a_ind = DG_MAT_IND(k, i, DG_NP_N1, DG_NP_N1);
        // int b_ind = fmaskL[j] * DG_NP_N1 + k;
        int b_ind = DG_MAT_IND(k, fmaskL_ind, DG_NP_N1, DG_NP_N1);
        op2L[op_ind] -= -0.5 * DL[a_ind] * sJ[0] * mmFL[b_ind];
      }
    }
  }

  // Do right face
  op2_in_kernel_gemm(false, false, DG_NP_N1, DG_NP_N1, DG_NP_N1, -0.5 * sJ[1], mmFR, DG_NP_N1, DR, DG_NP_N1, 1.0, op1R, DG_NP_N1);
  op2_in_kernel_gemm(true, false, DG_NP_N1, DG_NP_N1, DG_NP_N1, -0.5 * sJ[1], DR, DG_NP_N1, mmFR, DG_NP_N1, 1.0, op1R, DG_NP_N1);

  const DG_FP tmp_constR = 0.5 * sJ[1] * gtau;
  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int op_ind = i + j * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);
      op1R[op_ind] += tmp_constR * mmFR[op_ind];
    }
  }

  for(int i = 0; i < DG_NP_N1 * DG_NP_N1; i++) {
    op2R[i] = 0.0;
  }

  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NPF_N1; j++) {
      const int fmaskL_ind = *reverse ? fmaskL[DG_NPF_N1 - j - 1] : fmaskL[j];
      const int fmaskR_ind = fmaskR[j];
      // int op_ind = i + fmaskL_corrected[j] * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, fmaskL_ind, DG_NP_N1, DG_NP_N1);
      // int find = i + fmaskR[j] * DG_NP_N1;
      int find = DG_MAT_IND(i, fmaskR_ind, DG_NP_N1, DG_NP_N1);
      op2R[op_ind] -= 0.5 * gtau * sJ[1] * mmFR[find];
    }
  }

  for(int i = 0; i < DG_NPF_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      const int fmaskR_ind_i = fmaskR[i];
      // int op_ind = fmaskR[i] + j * DG_NP_N1;
      int op_ind = DG_MAT_IND(fmaskR_ind_i, j, DG_NP_N1, DG_NP_N1);
      for(int k = 0; k < DG_NPF_N1; k++) {
        const int fmaskR_ind_k = fmaskR[k];
        const int fmaskL_ind = *reverse ? fmaskL[DG_NPF_N1 - k - 1] : fmaskL[k];
        // int a_ind = fmaskR[i] + fmaskR[k] * DG_NP_N1;
        int a_ind = DG_MAT_IND(fmaskR_ind_i, fmaskR_ind_k, DG_NP_N1, DG_NP_N1);
        // int b_ind = j * DG_NP_N1 + fmaskL_corrected[k];
        int b_ind = DG_MAT_IND(fmaskL_ind, j, DG_NP_N1, DG_NP_N1);
        op2R[op_ind] -= 0.5 * sJ[1] * mmFR[a_ind] * -DL[b_ind];
      }
    }
  }

  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NPF_N1; j++) {
      const int fmaskL_ind = *reverse ? fmaskL[DG_NPF_N1 - j - 1] : fmaskL[j];
      // int op_ind = i + fmaskL_corrected[j] * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, fmaskL_ind, DG_NP_N1, DG_NP_N1);
      for(int k = 0; k < DG_NP_N1; k++) {
        const int fmaskR_ind = fmaskR[j];
        // int a_ind = i * DG_NP_N1 + k;
        int a_ind = DG_MAT_IND(k, i, DG_NP_N1, DG_NP_N1);
        // int b_ind = fmaskR[j] * DG_NP_N1 + k;
        int b_ind = DG_MAT_IND(k, fmaskR_ind, DG_NP_N1, DG_NP_N1);
        op2R[op_ind] -= -0.5 * DR[a_ind] * sJ[1] * mmFR[b_ind];
      }
    }
  }
}
