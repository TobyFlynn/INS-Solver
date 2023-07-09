inline void factor_poisson_coarse_matrix_2d_opbc(const int *faceNum,
                    const int *bc_type, const DG_FP *nx, const DG_FP *ny,
                    const DG_FP *fscale, const DG_FP *sJ, const DG_FP *geof,
                    const DG_FP *factor, DG_FP *op) {
  const DG_FP *dr_mat = dg_Dr_kernel;
  const DG_FP *ds_mat = dg_Ds_kernel;
  const DG_FP *mmF0_mat = dg_MM_F0_kernel;
  const DG_FP *mmF1_mat = dg_MM_F1_kernel;
  const DG_FP *mmF2_mat = dg_MM_F2_kernel;

  const DG_FP tau_order = 1.0; // (DG_FP) DG_ORDER;

  const DG_FP *mmF;
  if(*faceNum == 0)
    mmF = mmF0_mat;
  else if(*faceNum == 1)
    mmF = mmF1_mat;
  else
    mmF = mmF2_mat;

  const int find = *faceNum * DG_NPF_N1;
  const int *fmask  = FMASK;
  const int *fmaskB = &fmask[*faceNum * DG_NPF_N1];

  if(*bc_type == 0) {
    // Dirichlet
    DG_FP D[DG_NP_N1 * DG_NP_N1];
    const DG_FP r_fact = nx[0] * geof[RX_IND] + ny[0] * geof[RY_IND];
    const DG_FP s_fact = nx[0] * geof[SX_IND] + ny[0] * geof[SY_IND];
    for(int i = 0; i < DG_NP_N1; i++) {
      for(int j = 0; j < DG_NP_N1; j++) {
        // int ind = i + j * DG_NP_N1;
        int ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);
        D[ind] = r_fact * dr_mat[ind] + s_fact * ds_mat[ind];
        D[ind] *= factor[i];
      }
    }

    DG_FP gtau = 0.0;
    for(int i = 0; i < DG_NPF_N1; i++) {
      const int fmask_ind = fmaskB[i];
      gtau = fmax(gtau, (DG_FP)2.0 * (tau_order + 1) * (tau_order + 2) * *fscale * factor[fmask_ind]);
    }

    for(int i = 0; i < DG_NP_N1; i++) {
      for(int j = 0; j < DG_NPF_N1; j++) {
        const int fmask_ind = fmaskB[j];
        // int op_ind = i + j * DG_NP_N1;
        int op_ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);
        // int mm_ind = i + fmaskB[j] * DG_NP_N1;
        int mm_ind = DG_MAT_IND(i, fmask_ind, DG_NP_N1, DG_NP_N1);
        op[op_ind] = gtau * *sJ * mmF[mm_ind];
      }
    }

    for(int i = 0; i < DG_NP_N1; i++) {
      for(int j = 0; j < DG_NPF_N1; j++) {
        // int op_ind = i + j * DG_NP_N1;
        int op_ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NPF_N1);
        for(int k = 0; k < DG_NP_N1; k++) {
          const int fmask_ind = fmaskB[j];
          // int a_ind = i * DG_NP_N1 + k;
          int a_ind = DG_MAT_IND(k, i, DG_NP_N1, DG_NP_N1);
          // int b_ind = fmaskB[j] * DG_NP_N1 + k;
          int b_ind = DG_MAT_IND(k, fmask_ind, DG_NP_N1, DG_NP_N1);
          op[op_ind] += -D[a_ind] * *sJ * mmF[b_ind];
        }
      }
    }
  } else {
    // Neumann
    for(int i = 0; i < DG_NP_N1; i++) {
      for(int j = 0; j < DG_NPF_N1; j++) {
        const int fmask_ind = fmaskB[j];
        // int op_ind = i + j * DG_NP_N1;
        int op_ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NPF_N1);
        // int mm_ind = i + fmaskB[j] * DG_NP_N1;
        int mm_ind = DG_MAT_IND(i, fmask_ind, DG_NP_N1, DG_NP_N1);
        op[op_ind] = *sJ * mmF[mm_ind];
      }
    }
  }
}
