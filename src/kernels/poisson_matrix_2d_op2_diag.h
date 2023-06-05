inline void poisson_matrix_2d_op2_diag(const int *order, const int *faceNum,
                                       const bool *reversed, const DG_FP *nx,
                                       const DG_FP *ny, const DG_FP *fscale,
                                       const DG_FP *sJ, const DG_FP **geof,
                                       DG_FP *diagL, DG_FP *diagR) {
  const DG_FP *dr_mat = dg_Dr_kernel + (*order - 1) * DG_NP * DG_NP;
  const DG_FP *ds_mat = dg_Ds_kernel + (*order - 1) * DG_NP * DG_NP;

  const DG_FP *mmF0_mat = dg_MM_F0_kernel + (*order - 1) * DG_NP * DG_NP;
  const DG_FP *mmF1_mat = dg_MM_F1_kernel + (*order - 1) * DG_NP * DG_NP;
  const DG_FP *mmF2_mat = dg_MM_F2_kernel + (*order - 1) * DG_NP * DG_NP;

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

  const int dg_np  = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS + 1];
  const int findL = faceNum[0] * dg_npf;
  const int findR = faceNum[1] * dg_npf;
  const int *fmask  = FMASK;
  const int *fmaskL = &fmask[faceNum[0] * dg_npf];
  const int *fmaskR = &fmask[faceNum[1] * dg_npf];

  DG_FP DL[DG_NP * DG_NP], DR[DG_NP * DG_NP];
  const DG_FP r_fact_0 = nx[0] * geof[0][RX_IND] + ny[0] * geof[0][RY_IND];
  const DG_FP s_fact_0 = nx[0] * geof[0][SX_IND] + ny[0] * geof[0][SY_IND];
  const DG_FP r_fact_1 = nx[1] * geof[1][RX_IND] + ny[1] * geof[1][RY_IND];
  const DG_FP s_fact_1 = nx[1] * geof[1][SX_IND] + ny[1] * geof[1][SY_IND];
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int ind = i + j * DG_NP_N1;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      DL[ind] = r_fact_0 * dr_mat[ind] + s_fact_0 * ds_mat[ind];
      DR[ind] = r_fact_1 * dr_mat[ind] + s_fact_1 * ds_mat[ind];
    }
  }

  const DG_FP gtau = 2.0 * (*order + 1) * (*order + 2) * fmax(fscale[0], fscale[1]);

  for(int i = 0; i < dg_np; i++) {
    DG_FP tmp = 0.0;
    for(int k = 0; k < dg_np; k++) {
      int a_ind_t1 = DG_MAT_IND(i, k, dg_np, dg_np);
      int b_ind_t1 = DG_MAT_IND(k, i, dg_np, dg_np);
      int a_ind_t2 = DG_MAT_IND(k, i, dg_np, dg_np);
      int b_ind_t2 = DG_MAT_IND(k, i, dg_np, dg_np);
      tmp += mmFL[a_ind_t1] * DL[b_ind_t1] + DL[a_ind_t2] * mmFL[b_ind_t2];
    }
    int ind_t3 = DG_MAT_IND(i, i, dg_np, dg_np);
    diagL[i] += 0.5 * sJ[0] * (gtau * mmFL[ind_t3] - tmp);
  }

  for(int i = 0; i < dg_np; i++) {
    DG_FP tmp = 0.0;
    for(int k = 0; k < dg_np; k++) {
      int a_ind_t1 = DG_MAT_IND(i, k, dg_np, dg_np);
      int b_ind_t1 = DG_MAT_IND(k, i, dg_np, dg_np);
      int a_ind_t2 = DG_MAT_IND(k, i, dg_np, dg_np);
      int b_ind_t2 = DG_MAT_IND(k, i, dg_np, dg_np);
      tmp += mmFR[a_ind_t1] * DR[b_ind_t1] + DR[a_ind_t2] * mmFR[b_ind_t2];
    }
    int ind_t3 = DG_MAT_IND(i, i, dg_np, dg_np);
    diagR[i] += 0.5 * sJ[1] * (gtau * mmFR[ind_t3] - tmp);
  }
}
