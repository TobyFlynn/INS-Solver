inline void poisson_coarse_matrix_2d_op2(const int *faceNum, const bool *reversed,
                                         const DG_FP *nx, const DG_FP *ny,
                                         const DG_FP *fscale, const DG_FP *sJ,
                                         const DG_FP **geof, DG_FP *op1L, DG_FP *op1R,
                                         DG_FP *op2L, DG_FP *op2R) {
  const DG_FP *dr_mat = dg_Dr_kernel;
  const DG_FP *ds_mat = dg_Ds_kernel;

  const DG_FP *mmF0_mat = dg_MM_F0_kernel;
  const DG_FP *mmF1_mat = dg_MM_F1_kernel;
  const DG_FP *mmF2_mat = dg_MM_F2_kernel;

  const bool rev = *reversed;

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

  const int findL = faceNum[0] * DG_NPF_N1;
  const int findR = faceNum[1] * DG_NPF_N1;
  const int *fmask  = FMASK;
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF_N1];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF_N1];

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
      DR[ind] = r_fact_1 * dr_mat[ind] + s_fact_1 * ds_mat[ind];
    }
  }

  const DG_FP gtau = 2.0 * (1 + 1) * (1 + 2) * fmax(fscale[0], fscale[1]);

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
      // int op_ind = i + fmaskR_corrected[j] * DG_NP_N1;
      const int fmaskR_ind = rev ? fmaskR[DG_NPF_N1 - j - 1] : fmaskR[j];
      int op_ind = DG_MAT_IND(i, fmaskR_ind, DG_NP_N1, DG_NP_N1);
      // int find = i + fmaskL[j] * DG_NP_N1;
      int find = DG_MAT_IND(i, fmaskL[j], DG_NP_N1, DG_NP_N1);
      op2L[op_ind] -= gtau * mmFL[find];
    }
  }

  for(int i = 0; i < DG_NPF_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int op_ind = fmaskL[i] + j * DG_NP_N1;
      int op_ind = DG_MAT_IND(fmaskL[i], j, DG_NP_N1, DG_NP_N1);
      for(int k = 0; k < DG_NPF_N1; k++) {
        // int a_ind = fmaskL[i] + fmaskL[k] * DG_NP_N1;
        int a_ind = DG_MAT_IND(fmaskL[i], fmaskL[k], DG_NP_N1, DG_NP_N1);
        // int b_ind = j * DG_NP_N1 + fmaskR_corrected[k];
        const int fmaskR_ind = rev ? fmaskR[DG_NPF_N1 - k - 1] : fmaskR[k];
        int b_ind = DG_MAT_IND(fmaskR_ind, j, DG_NP_N1, DG_NP_N1);
        op2L[op_ind] += mmFL[a_ind] * DR[b_ind];
      }
    }
  }

  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NPF_N1; j++) {
      const int fmaskR_ind = rev ? fmaskR[DG_NPF_N1 - j - 1] : fmaskR[j];
      // int op_ind = i + fmaskR_corrected[j] * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, fmaskR_ind, DG_NP_N1, DG_NP_N1);
      for(int k = 0; k < DG_NP_N1; k++) {
        // int a_ind = i * DG_NP_N1 + k;
        int a_ind = DG_MAT_IND(k, i, DG_NP_N1, DG_NP_N1);
        // int b_ind = fmaskL[j] * DG_NP_N1 + k;
        int b_ind = DG_MAT_IND(k, fmaskL[j], DG_NP_N1, DG_NP_N1);
        op2L[op_ind] += DL[a_ind] * mmFL[b_ind];
      }
    }
  }

  for(int i = 0; i < DG_NP_N1 * DG_NP_N1; i++) {
    op2L[i] *= 0.5 * sJ[0];
  }

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
      const int fmaskL_ind = rev ? fmaskL[DG_NPF_N1 - j - 1] : fmaskL[j];
      // int op_ind = i + fmaskL_corrected[j] * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, fmaskL_ind, DG_NP_N1, DG_NP_N1);
      // int find = i + fmaskR[j] * DG_NP_N1;
      int find = DG_MAT_IND(i, fmaskR[j], DG_NP_N1, DG_NP_N1);
      op2R[op_ind] -= gtau * mmFR[find];
    }
  }

  for(int i = 0; i < DG_NPF_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int op_ind = fmaskR[i] + j * DG_NP_N1;
      int op_ind = DG_MAT_IND(fmaskR[i], j, DG_NP_N1, DG_NP_N1);
      for(int k = 0; k < DG_NPF_N1; k++) {
        // int a_ind = fmaskR[i] + fmaskR[k] * DG_NP_N1;
        int a_ind = DG_MAT_IND(fmaskR[i], fmaskR[k], DG_NP_N1, DG_NP_N1);
        // int b_ind = j * DG_NP_N1 + fmaskL_corrected[k];
        const int fmaskL_ind = rev ? fmaskL[DG_NPF_N1 - k - 1] : fmaskL[k];
        int b_ind = DG_MAT_IND(fmaskL_ind, j, DG_NP_N1, DG_NP_N1);
        op2R[op_ind] += mmFR[a_ind] * DL[b_ind];
      }
    }
  }

  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NPF_N1; j++) {
      const int fmaskL_ind = rev ? fmaskL[DG_NPF_N1 - j - 1] : fmaskL[j];
      // int op_ind = i + fmaskL_corrected[j] * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, fmaskL_ind, DG_NP_N1, DG_NP_N1);
      for(int k = 0; k < DG_NP_N1; k++) {
        // int a_ind = i * DG_NP_N1 + k;
        int a_ind = DG_MAT_IND(k, i, DG_NP_N1, DG_NP_N1);
        // int b_ind = fmaskR[j] * DG_NP_N1 + k;
        int b_ind = DG_MAT_IND(k, fmaskR[j], DG_NP_N1, DG_NP_N1);
        op2R[op_ind] += DR[a_ind] * mmFR[b_ind];
      }
    }
  }

  for(int i = 0; i < DG_NP_N1 * DG_NP_N1; i++) {
    op2R[i] *= 0.5 * sJ[1];
  }
}
