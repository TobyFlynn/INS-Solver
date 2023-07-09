inline void factor_poisson_matrix_2d_op2_partial_diag(const int *order,
                  const int *faceNum, const bool *reverse, const DG_FP *nx,
                  const DG_FP *ny, const DG_FP *fscale, const DG_FP *sJ,
                  const DG_FP **geof, const DG_FP **factor, DG_FP *diagL,
                  DG_FP *diagR) {
  const int p = *order;
  const DG_FP *dr_mat = &dg_Dr_kernel[(p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &dg_Ds_kernel[(p - 1) * DG_NP * DG_NP];
  const DG_FP *mmF0_mat = &dg_MM_F0_kernel[(p - 1) * DG_NP * DG_NP];
  const DG_FP *mmF1_mat = &dg_MM_F1_kernel[(p - 1) * DG_NP * DG_NP];
  const DG_FP *mmF2_mat = &dg_MM_F2_kernel[(p - 1) * DG_NP * DG_NP];
  const int dg_np  = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

  const DG_FP tau_order = (DG_FP) p; // (DG_FP) DG_ORDER;

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
  const int findL = faceNumL * dg_npf;
  const int findR = faceNumR * dg_npf;
  const int *fmask  = &FMASK[(p - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskL = &fmask[faceNumL * dg_npf];
  const int *fmaskR = &fmask[faceNumR * dg_npf];

  DG_FP D[DG_NP * DG_NP];
  const DG_FP r_fact_0 = nx[0] * geof[0][RX_IND] + ny[0] * geof[0][RY_IND];
  const DG_FP s_fact_0 = nx[0] * geof[0][SX_IND] + ny[0] * geof[0][SY_IND];
  for(int j = 0; j < dg_np; j++) {
    for(int i = 0; i < dg_np; i++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      D[ind] = r_fact_0 * dr_mat[ind] + s_fact_0 * ds_mat[ind];
      D[ind] *= factor[0][i];
    }
  }

  DG_FP gtau = 0.0;
  for(int i = 0; i < dg_npf; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = *reverse ? fmaskR[dg_npf - i - 1] : fmaskR[i];
    DG_FP tmp = 2.0 * (tau_order + 1) * (tau_order + 2) * fmax(fscale[0] * factor[0][fmaskL_ind], fscale[1] * factor[1][fmaskR_ind]);
    gtau = fmax(gtau, tmp);
  }

  for(int i = 0; i < dg_np; i++) {
    DG_FP tmp = 0.0;
    for(int k = 0; k < dg_np; k++) {
      int a_ind_t1 = DG_MAT_IND(i, k, dg_np, dg_np);
      int b_ind_t1 = DG_MAT_IND(k, i, dg_np, dg_np);
      int a_ind_t2 = DG_MAT_IND(k, i, dg_np, dg_np);
      int b_ind_t2 = DG_MAT_IND(k, i, dg_np, dg_np);
      tmp += mmFL[a_ind_t1] * D[b_ind_t1] + D[a_ind_t2] * mmFL[b_ind_t2];
    }
    int ind_t3 = DG_MAT_IND(i, i,  dg_np,  dg_np);
    diagL[i] += 0.5 * sJ[0] * (gtau * mmFL[ind_t3] - tmp);
  }

  const DG_FP r_fact_1 = nx[1] * geof[1][RX_IND] + ny[1] * geof[1][RY_IND];
  const DG_FP s_fact_1 = nx[1] * geof[1][SX_IND] + ny[1] * geof[1][SY_IND];
  for(int j = 0; j < dg_np; j++) {
    for(int i = 0; i < dg_np; i++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      D[ind] = r_fact_1 * dr_mat[ind] + s_fact_1 * ds_mat[ind];
      D[ind] *= factor[1][i];
    }
  }

  for(int i = 0; i < dg_np; i++) {
    DG_FP tmp = 0.0;
    for(int k = 0; k < dg_np; k++) {
      int a_ind_t1 = DG_MAT_IND(i, k, dg_np, dg_np);
      int b_ind_t1 = DG_MAT_IND(k, i, dg_np, dg_np);
      int a_ind_t2 = DG_MAT_IND(k, i, dg_np, dg_np);
      int b_ind_t2 = DG_MAT_IND(k, i, dg_np, dg_np);
      tmp += mmFR[a_ind_t1] * D[b_ind_t1] + D[a_ind_t2] * mmFR[b_ind_t2];
    }
    int ind_t3 = DG_MAT_IND(i, i, dg_np, dg_np);
    diagR[i] += 0.5 * sJ[1] * (gtau * mmFR[ind_t3] - tmp);
  }
}
