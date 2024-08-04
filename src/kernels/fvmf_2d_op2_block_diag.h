inline void fvmf_2d_op2_block_diag(const int *faceNum, const bool *reverse, 
                  const DG_FP *nx, const DG_FP *ny, const DG_FP *fscale, 
                  const DG_FP *sJ, const DG_FP **geof, const DG_FP **factor, 
                  DG_FP *op1L, DG_FP *op1R) {
  const DG_FP *dr_mat = &dg_Dr_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &dg_Ds_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF0_mat = &dg_MM_F0_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF1_mat = &dg_MM_F1_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF2_mat = &dg_MM_F2_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];

  const DG_FP tau_order = (DG_FP) DG_ORDER;

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
  const int findL = faceNumL * DG_NPF;
  const int findR = faceNumR * DG_NPF;
  const int *fmask  = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskL = &fmask[faceNumL * DG_NPF];
  const int *fmaskR = &fmask[faceNumR * DG_NPF];

  DG_FP D[DG_NP * DG_NP];
  const DG_FP r_fact_0 = nx[0] * geof[0][RX_IND] + ny[0] * geof[0][RY_IND];
  const DG_FP s_fact_0 = nx[0] * geof[0][SX_IND] + ny[0] * geof[0][SY_IND];
  for(int j = 0; j < DG_NP; j++) {
    for(int i = 0; i < DG_NP; i++) {
      // int ind = i + j * DG_NP;
      int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      D[ind] = r_fact_0 * dr_mat[ind] + s_fact_0 * ds_mat[ind];
      D[ind] *= factor[0][i];
    }
  }

  DG_FP gtau = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = *reverse ? fmaskR[DG_NPF - i - 1] : fmaskR[i];
    DG_FP tmp = 2.0 * (tau_order + 1) * (tau_order + 2) * fmax(fscale[0] * factor[0][fmaskL_ind], fscale[1] * factor[1][fmaskR_ind]);
    gtau = fmax(gtau, tmp);
  }

  op2_in_kernel_gemm(false, false, DG_NP, DG_NP, DG_NP, -0.5 * sJ[0], mmFL, DG_NP, D, DG_NP, 1.0, op1L, DG_NP);
  op2_in_kernel_gemm(true, false, DG_NP, DG_NP, DG_NP, -0.5 * sJ[0], D, DG_NP, mmFL, DG_NP, 1.0, op1L, DG_NP);

  const DG_FP tmp_constL = 0.5 * sJ[0] * gtau;
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      // int op_ind = i + j * DG_NP;
      int op_ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      op1L[op_ind] += tmp_constL * mmFL[op_ind];
    }
  }

  const DG_FP r_fact_1 = nx[1] * geof[1][RX_IND] + ny[1] * geof[1][RY_IND];
  const DG_FP s_fact_1 = nx[1] * geof[1][SX_IND] + ny[1] * geof[1][SY_IND];
  for(int j = 0; j < DG_NP; j++) {
    for(int i = 0; i < DG_NP; i++) {
      // int ind = i + j * DG_NP;
      int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      D[ind] = r_fact_1 * dr_mat[ind] + s_fact_1 * ds_mat[ind];
      D[ind] *= factor[1][i];
    }
  }

  op2_in_kernel_gemm(false, false, DG_NP, DG_NP, DG_NP, -0.5 * sJ[1], mmFR, DG_NP, D, DG_NP, 1.0, op1R, DG_NP);
  op2_in_kernel_gemm(true, false, DG_NP, DG_NP, DG_NP, -0.5 * sJ[1], D, DG_NP, mmFR, DG_NP, 1.0, op1R, DG_NP);

  const DG_FP tmp_constR = 0.5 * sJ[1] * gtau;
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      // int op_ind = i + j * DG_NP;
      int op_ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      op1R[op_ind] += tmp_constR * mmFR[op_ind];
    }
  }
}
