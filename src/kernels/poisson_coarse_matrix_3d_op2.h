inline void poisson_coarse_matrix_3d_op2(const int *faceNum, const int *fmaskL_corrected,
                                         const int *fmaskR_corrected, const DG_FP *nx,
                                         const DG_FP *ny, const DG_FP *nz,
                                         const DG_FP *fscale, const DG_FP *sJ,
                                         const DG_FP **rx, const DG_FP **sx,
                                         const DG_FP **tx, const DG_FP **ry,
                                         const DG_FP **sy, const DG_FP **ty,
                                         const DG_FP **rz, const DG_FP **sz,
                                         const DG_FP **tz, DG_FP *op1L, DG_FP *op1R,
                                         DG_FP *op2L, DG_FP *op2R) {
  const DG_FP *dr_mat = dg_Dr_kernel;
  const DG_FP *ds_mat = dg_Ds_kernel;
  const DG_FP *dt_mat = dg_Dt_kernel;
  const DG_FP *mmF0_mat = dg_MM_F0_kernel;
  const DG_FP *mmF1_mat = dg_MM_F1_kernel;
  const DG_FP *mmF2_mat = dg_MM_F2_kernel;
  const DG_FP *mmF3_mat = dg_MM_F3_kernel;

  const DG_FP *mmFL, *mmFR;
  if(faceNum[0] == 0)
    mmFL = mmF0_mat;
  else if(faceNum[0] == 1)
    mmFL = mmF1_mat;
  else if(faceNum[0] == 2)
    mmFL = mmF2_mat;
  else
    mmFL = mmF3_mat;

  if(faceNum[1] == 0)
    mmFR = mmF0_mat;
  else if(faceNum[1] == 1)
    mmFR = mmF1_mat;
  else if(faceNum[1] == 2)
    mmFR = mmF2_mat;
  else
    mmFR = mmF3_mat;

  const int faceNumL = faceNum[0];
  const int faceNumR = faceNum[1];
  const int findL = faceNumL * DG_NPF_N1;
  const int findR = faceNumR * DG_NPF_N1;
  const int *fmask  = FMASK;
  const int *fmaskL = &fmask[faceNumL * DG_NPF_N1];
  const int *fmaskR = &fmask[faceNumR * DG_NPF_N1];

  DG_FP DL[DG_NP_N1 * DG_NP_N1], DR[DG_NP_N1 * DG_NP_N1];
  const DG_FP r_fact_0 = nx[0] * rx[0][0] + ny[0] * ry[0][0] + nz[0] * rz[0][0];
  const DG_FP s_fact_0 = nx[0] * sx[0][0] + ny[0] * sy[0][0] + nz[0] * sz[0][0];
  const DG_FP t_fact_0 = nx[0] * tx[0][0] + ny[0] * ty[0][0] + nz[0] * tz[0][0];
  const DG_FP r_fact_1 = nx[1] * rx[1][0] + ny[1] * ry[1][0] + nz[1] * rz[1][0];
  const DG_FP s_fact_1 = nx[1] * sx[1][0] + ny[1] * sy[1][0] + nz[1] * sz[1][0];
  const DG_FP t_fact_1 = nx[1] * tx[1][0] + ny[1] * ty[1][0] + nz[1] * tz[1][0];
  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int ind = i + j * DG_NP_N1;
      int ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);
      DL[ind] = r_fact_0 * dr_mat[ind] + s_fact_0 * ds_mat[ind] + t_fact_0 * dt_mat[ind];
      DR[ind] = r_fact_1 * dr_mat[ind] + s_fact_1 * ds_mat[ind] + t_fact_1 * dt_mat[ind];
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
      const int fmaskL_ind = fmaskL[j];
      const int fmaskR_ind = fmaskR_corrected[j];
      // int op_ind = i + fmaskR_corrected[j] * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, fmaskR_ind, DG_NP_N1, DG_NP_N1);
      // int find = i + fmaskL[j] * DG_NP_N1;
      int find = DG_MAT_IND(i, fmaskL_ind, DG_NP_N1, DG_NP_N1);
      op2L[op_ind] -= gtau * mmFL[find];
    }
  }

  for(int i = 0; i < DG_NPF_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      const int fmaskL_ind_i = fmaskL[i];
      // int op_ind = fmaskL[i] + j * DG_NP_N1;
      int op_ind = DG_MAT_IND(fmaskL_ind_i, j, DG_NP_N1, DG_NP_N1);
      for(int k = 0; k < DG_NPF_N1; k++) {
        const int fmaskL_ind_k = fmaskL[k];
        const int fmaskR_ind_k = fmaskR_corrected[k];
        // int a_ind = fmaskL[i] + fmaskL[k] * DG_NP_N1;
        int a_ind = DG_MAT_IND(fmaskL_ind_i, fmaskL_ind_k, DG_NP_N1, DG_NP_N1);
        // int b_ind = j * DG_NP_N1 + fmaskR_corrected[k];
        int b_ind = DG_MAT_IND(fmaskR_ind_k, j, DG_NP_N1, DG_NP_N1);
        op2L[op_ind] += mmFL[a_ind] * DR[b_ind];
      }
    }
  }

  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NPF_N1; j++) {
      const int fmaskR_ind = fmaskR_corrected[j];
      // int op_ind = i + fmaskR_corrected[j] * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, fmaskR_ind, DG_NP_N1, DG_NP_N1);
      for(int k = 0; k < DG_NP_N1; k++) {
        const int fmaskL_ind = fmaskL[j];
        // int a_ind = i * DG_NP_N1 + k;
        int a_ind = DG_MAT_IND(k, i, DG_NP_N1, DG_NP_N1);
        // int b_ind = fmaskL[j] * DG_NP_N1 + k;
        int b_ind = DG_MAT_IND(k, fmaskL_ind, DG_NP_N1, DG_NP_N1);
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
      const int fmaskL_ind = fmaskL_corrected[j];
      const int fmaskR_ind = fmaskR[j];
      // int op_ind = i + fmaskL_corrected[j] * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, fmaskL_ind, DG_NP_N1, DG_NP_N1);
      // int find = i + fmaskR[j] * DG_NP_N1;
      int find = DG_MAT_IND(i, fmaskR_ind, DG_NP_N1, DG_NP_N1);
      op2R[op_ind] -= gtau * mmFR[find];
    }
  }

  for(int i = 0; i < DG_NPF_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      const int fmaskR_ind_i = fmaskR[i];
      // int op_ind = fmaskR[i] + j * DG_NP_N1;
      int op_ind = DG_MAT_IND(fmaskR_ind_i, j, DG_NP_N1, DG_NP_N1);
      for(int k = 0; k < DG_NPF_N1; k++) {
        const int fmaskR_ind_k = fmaskR[k];
        const int fmaskL_ind = fmaskL_corrected[k];
        // int a_ind = fmaskR[i] + fmaskR[k] * DG_NP_N1;
        int a_ind = DG_MAT_IND(fmaskR_ind_i, fmaskR_ind_k, DG_NP_N1, DG_NP_N1);
        // int b_ind = j * DG_NP_N1 + fmaskL_corrected[k];
        int b_ind = DG_MAT_IND(fmaskL_ind, j, DG_NP_N1, DG_NP_N1);
        op2R[op_ind] += mmFR[a_ind] * DL[b_ind];
      }
    }
  }

  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NPF_N1; j++) {
      const int fmaskL_ind = fmaskL_corrected[j];
      // int op_ind = i + fmaskL_corrected[j] * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, fmaskL_ind, DG_NP_N1, DG_NP_N1);
      for(int k = 0; k < DG_NP_N1; k++) {
        const int fmaskR_ind = fmaskR[j];
        // int a_ind = i * DG_NP_N1 + k;
        int a_ind = DG_MAT_IND(k, i, DG_NP_N1, DG_NP_N1);
        // int b_ind = fmaskR[j] * DG_NP_N1 + k;
        int b_ind = DG_MAT_IND(k, fmaskR_ind, DG_NP_N1, DG_NP_N1);
        op2R[op_ind] += DR[a_ind] * mmFR[b_ind];
      }
    }
  }

  for(int i = 0; i < DG_NP_N1 * DG_NP_N1; i++) {
    op2R[i] *= 0.5 * sJ[1];
  }
}
