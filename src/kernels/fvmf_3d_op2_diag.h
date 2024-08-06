inline void fvmf_3d_op2_diag(const int *faceNum, const int *fmaskL_corrected, const int *fmaskR_corrected,
                  const DG_FP *nx, const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale, const DG_FP *sJ, 
                  const DG_FP **geof, const DG_FP **factor, DG_FP *diagL, DG_FP *diagR) {
  const DG_FP *dr_mat = &dg_Dr_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &dg_Ds_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dg_Dt_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF0_mat = &dg_MM_F0_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF1_mat = &dg_MM_F1_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF2_mat = &dg_MM_F2_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF3_mat = &dg_MM_F3_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];

  const DG_FP tau_order = (DG_FP) DG_ORDER;

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
  const int findL = faceNumL * DG_NPF;
  const int findR = faceNumR * DG_NPF;
  const int *fmask  = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskL = &fmask[faceNumL * DG_NPF];
  const int *fmaskR = &fmask[faceNumR * DG_NPF];

  DG_FP D[DG_NP * DG_NP];
  const DG_FP r_fact_0 = nx[0] * geof[0][RX_IND] + ny[0] * geof[0][RY_IND] + nz[0] * geof[0][RZ_IND];
  const DG_FP s_fact_0 = nx[0] * geof[0][SX_IND] + ny[0] * geof[0][SY_IND] + nz[0] * geof[0][SZ_IND];
  const DG_FP t_fact_0 = nx[0] * geof[0][TX_IND] + ny[0] * geof[0][TY_IND] + nz[0] * geof[0][TZ_IND];
  for(int j = 0; j < DG_NP; j++) {
    for(int i = 0; i < DG_NP; i++) {
      // int ind = i + j * DG_NP;
      int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      D[ind] = r_fact_0 * dr_mat[ind] + s_fact_0 * ds_mat[ind] + t_fact_0 * dt_mat[ind];
      D[ind] *= factor[0][i];
    }
  }

  DG_FP gtau = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];
    DG_FP tmp = 2.0 * (tau_order + 1) * (tau_order + 2) * fmax(fscale[0] * factor[0][fmaskL_ind], fscale[1] * factor[1][fmaskR_ind]);
    gtau = fmax(gtau, tmp);
  }

  for(int i = 0; i < DG_NP; i++) {
    DG_FP tmp = 0.0;
    for(int k = 0; k < DG_NP; k++) {
      int a_ind_t1 = DG_MAT_IND(i, k, DG_NP, DG_NP);
      int b_ind_t1 = DG_MAT_IND(k, i, DG_NP, DG_NP);
      int a_ind_t2 = DG_MAT_IND(k, i, DG_NP, DG_NP);
      int b_ind_t2 = DG_MAT_IND(k, i, DG_NP, DG_NP);
      tmp += mmFL[a_ind_t1] * D[b_ind_t1] + D[a_ind_t2] * mmFL[b_ind_t2];
    }
    int ind_t3 = DG_MAT_IND(i, i,  DG_NP,  DG_NP);
    diagL[i] += 0.5 * sJ[0] * (gtau * mmFL[ind_t3] - tmp);
  }

  const DG_FP r_fact_1 = nx[1] * geof[1][RX_IND] + ny[1] * geof[1][RY_IND] + nz[1] * geof[1][RZ_IND];
  const DG_FP s_fact_1 = nx[1] * geof[1][SX_IND] + ny[1] * geof[1][SY_IND] + nz[1] * geof[1][SZ_IND];
  const DG_FP t_fact_1 = nx[1] * geof[1][TX_IND] + ny[1] * geof[1][TY_IND] + nz[1] * geof[1][TZ_IND];
  for(int j = 0; j < DG_NP; j++) {
    for(int i = 0; i < DG_NP; i++) {
      // int ind = i + j * DG_NP;
      int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      D[ind] = r_fact_1 * dr_mat[ind] + s_fact_1 * ds_mat[ind] + t_fact_1 * dt_mat[ind];
      D[ind] *= factor[1][i];
    }
  }

  for(int i = 0; i < DG_NP; i++) {
    DG_FP tmp = 0.0;
    for(int k = 0; k < DG_NP; k++) {
      int a_ind_t1 = DG_MAT_IND(i, k, DG_NP, DG_NP);
      int b_ind_t1 = DG_MAT_IND(k, i, DG_NP, DG_NP);
      int a_ind_t2 = DG_MAT_IND(k, i, DG_NP, DG_NP);
      int b_ind_t2 = DG_MAT_IND(k, i, DG_NP, DG_NP);
      tmp += mmFR[a_ind_t1] * D[b_ind_t1] + D[a_ind_t2] * mmFR[b_ind_t2];
    }
    int ind_t3 = DG_MAT_IND(i, i, DG_NP, DG_NP);
    diagR[i] += 0.5 * sJ[1] * (gtau * mmFR[ind_t3] - tmp);
  }
}
