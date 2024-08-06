inline void fvmf_3d_bop_block_diag(const int *faceNum, const int *u_bc_type, 
                  const int *v_bc_type, const int *w_bc_type, const DG_FP *nx,
                  const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale, 
                  const DG_FP *sJ, const DG_FP *geof, const DG_FP *factor,
                  DG_FP *u_op1, DG_FP *v_op1, DG_FP *w_op1) {
  const DG_FP tau_order = (DG_FP) DG_ORDER;
  const DG_FP *dr_mat = &dg_Dr_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &dg_Ds_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dg_Dt_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF0_mat = &dg_MM_F0_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF1_mat = &dg_MM_F1_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF2_mat = &dg_MM_F2_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF3_mat = &dg_MM_F3_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF;
  if(*faceNum == 0)
    mmF = mmF0_mat;
  else if(*faceNum == 1)
    mmF = mmF1_mat;
  else if(*faceNum == 2)
    mmF = mmF2_mat;
  else
    mmF = mmF3_mat;
  const int find = *faceNum * DG_NPF;
  const int *fmask  = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * DG_NPF];

  if(*u_bc_type == 2 && *v_bc_type == 2 && *w_bc_type == 2) {
    DG_FP D[DG_NP * DG_NP];
    const DG_FP r_fact = nx[0] * geof[RX_IND] + ny[0] * geof[RY_IND] + nz[0] * geof[RZ_IND];
    const DG_FP s_fact = nx[0] * geof[SX_IND] + ny[0] * geof[SY_IND] + nz[0] * geof[SZ_IND];
    const DG_FP t_fact = nx[0] * geof[TX_IND] + ny[0] * geof[TY_IND] + nz[0] * geof[TZ_IND];
    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int ind = i + j * DG_NP;
        int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
        D[ind] = r_fact * dr_mat[ind] + s_fact * ds_mat[ind] + t_fact * dt_mat[ind];
        D[ind] *= factor[i];
      }
    }

    DG_FP gtau = 0.0;
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      gtau = fmax(gtau, (DG_FP)2.0 * (tau_order + 1) * (tau_order + 2) * *fscale * factor[fmask_ind]);
    }

    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int op_ind = i + j * DG_NP;
        int op_ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
        DG_FP tmp = 0.0;
        for(int k = 0; k < DG_NP; k++) {
          // int a_ind0 = i + k * DG_NP;
          int a_ind0 = DG_MAT_IND(i, k, DG_NP, DG_NP);
          // int a_ind1 = i * DG_NP + k;
          int a_ind1 = DG_MAT_IND(k, i, DG_NP, DG_NP);
          // int b_ind  = j * DG_NP + k;
          int b_ind  = DG_MAT_IND(k, j, DG_NP, DG_NP);
          tmp += mmF[a_ind0] * D[b_ind] + D[a_ind1] * mmF[b_ind];
        }
        u_op1[op_ind] += *sJ * *nx * *nx * (gtau * mmF[op_ind] - tmp);
        v_op1[op_ind] += *sJ * *ny * *ny * (gtau * mmF[op_ind] - tmp);
        w_op1[op_ind] += *sJ * *nz * *nz * (gtau * mmF[op_ind] - tmp);
      }
    }
    return;
  }
  
  if(*u_bc_type == 1) {
    // Do nothing for Neumann boundary conditions
  } else if(*u_bc_type == 0) {
    DG_FP D[DG_NP * DG_NP];
    const DG_FP r_fact = nx[0] * geof[RX_IND] + ny[0] * geof[RY_IND] + nz[0] * geof[RZ_IND];
    const DG_FP s_fact = nx[0] * geof[SX_IND] + ny[0] * geof[SY_IND] + nz[0] * geof[SZ_IND];
    const DG_FP t_fact = nx[0] * geof[TX_IND] + ny[0] * geof[TY_IND] + nz[0] * geof[TZ_IND];
    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int ind = i + j * DG_NP;
        int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
        D[ind] = r_fact * dr_mat[ind] + s_fact * ds_mat[ind] + t_fact * dt_mat[ind];
        D[ind] *= factor[i];
      }
    }

    DG_FP gtau = 0.0;
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      gtau = fmax(gtau, (DG_FP)2.0 * (tau_order + 1) * (tau_order + 2) * *fscale * factor[fmask_ind]);
    }

    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int op_ind = i + j * DG_NP;
        int op_ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
        DG_FP tmp = 0.0;
        for(int k = 0; k < DG_NP; k++) {
          // int a_ind0 = i + k * DG_NP;
          int a_ind0 = DG_MAT_IND(i, k, DG_NP, DG_NP);
          // int a_ind1 = i * DG_NP + k;
          int a_ind1 = DG_MAT_IND(k, i, DG_NP, DG_NP);
          // int b_ind  = j * DG_NP + k;
          int b_ind  = DG_MAT_IND(k, j, DG_NP, DG_NP);
          tmp += -*sJ * mmF[a_ind0] * D[b_ind] - D[a_ind1] * *sJ * mmF[b_ind];
        }
        u_op1[op_ind] += gtau * *sJ * mmF[op_ind] + tmp;
      }
    }
  }

  if(*v_bc_type == 1) {
    // Do nothing for Neumann boundary conditions
  } else if(*v_bc_type == 0) {
    DG_FP D[DG_NP * DG_NP];
    const DG_FP r_fact = nx[0] * geof[RX_IND] + ny[0] * geof[RY_IND] + nz[0] * geof[RZ_IND];
    const DG_FP s_fact = nx[0] * geof[SX_IND] + ny[0] * geof[SY_IND] + nz[0] * geof[SZ_IND];
    const DG_FP t_fact = nx[0] * geof[TX_IND] + ny[0] * geof[TY_IND] + nz[0] * geof[TZ_IND];
    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int ind = i + j * DG_NP;
        int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
        D[ind] = r_fact * dr_mat[ind] + s_fact * ds_mat[ind] + t_fact * dt_mat[ind];
        D[ind] *= factor[i];
      }
    }

    DG_FP gtau = 0.0;
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      gtau = fmax(gtau, (DG_FP)2.0 * (tau_order + 1) * (tau_order + 2) * *fscale * factor[fmask_ind]);
    }

    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int op_ind = i + j * DG_NP;
        int op_ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
        DG_FP tmp = 0.0;
        for(int k = 0; k < DG_NP; k++) {
          // int a_ind0 = i + k * DG_NP;
          int a_ind0 = DG_MAT_IND(i, k, DG_NP, DG_NP);
          // int a_ind1 = i * DG_NP + k;
          int a_ind1 = DG_MAT_IND(k, i, DG_NP, DG_NP);
          // int b_ind  = j * DG_NP + k;
          int b_ind  = DG_MAT_IND(k, j, DG_NP, DG_NP);
          tmp += -*sJ * mmF[a_ind0] * D[b_ind] - D[a_ind1] * *sJ * mmF[b_ind];
        }
        v_op1[op_ind] += gtau * *sJ * mmF[op_ind] + tmp;
      }
    }
  }

  if(*w_bc_type == 1) {
    // Do nothing for Neumann boundary conditions
  } else if(*w_bc_type == 0) {
    DG_FP D[DG_NP * DG_NP];
    const DG_FP r_fact = nx[0] * geof[RX_IND] + ny[0] * geof[RY_IND] + nz[0] * geof[RZ_IND];
    const DG_FP s_fact = nx[0] * geof[SX_IND] + ny[0] * geof[SY_IND] + nz[0] * geof[SZ_IND];
    const DG_FP t_fact = nx[0] * geof[TX_IND] + ny[0] * geof[TY_IND] + nz[0] * geof[TZ_IND];
    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int ind = i + j * DG_NP;
        int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
        D[ind] = r_fact * dr_mat[ind] + s_fact * ds_mat[ind] + t_fact * dt_mat[ind];
        D[ind] *= factor[i];
      }
    }

    DG_FP gtau = 0.0;
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      gtau = fmax(gtau, (DG_FP)2.0 * (tau_order + 1) * (tau_order + 2) * *fscale * factor[fmask_ind]);
    }

    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int op_ind = i + j * DG_NP;
        int op_ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
        DG_FP tmp = 0.0;
        for(int k = 0; k < DG_NP; k++) {
          // int a_ind0 = i + k * DG_NP;
          int a_ind0 = DG_MAT_IND(i, k, DG_NP, DG_NP);
          // int a_ind1 = i * DG_NP + k;
          int a_ind1 = DG_MAT_IND(k, i, DG_NP, DG_NP);
          // int b_ind  = j * DG_NP + k;
          int b_ind  = DG_MAT_IND(k, j, DG_NP, DG_NP);
          tmp += -*sJ * mmF[a_ind0] * D[b_ind] - D[a_ind1] * *sJ * mmF[b_ind];
        }
        w_op1[op_ind] += gtau * *sJ * mmF[op_ind] + tmp;
      }
    }
  }
}
