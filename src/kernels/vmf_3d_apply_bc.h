inline void vmf_3d_apply_bc(const int *faceNum, const int *u_bc_type, const int *v_bc_type, 
                  const int *w_bc_type, const DG_FP *nx, const DG_FP *ny, const DG_FP *nz, 
                  const DG_FP *fscale, const DG_FP *sJ, const DG_FP *geof, const DG_FP *u_bc, 
                  const DG_FP *v_bc, const DG_FP *w_bc, DG_FP *u_rhs, DG_FP *v_rhs, DG_FP *w_rhs) {
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

  // u bcs
  if(*u_bc_type == 0) {
    // Dirichlet
    DG_FP D[DG_NP * DG_NP];
    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int ind = i + j * DG_NP;
        int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);

        D[ind] = *nx * (geof[RX_IND] * dr_mat[ind] + geof[SX_IND] * ds_mat[ind] + geof[TX_IND] * dt_mat[ind]);
        D[ind] += *ny * (geof[RY_IND] * dr_mat[ind] + geof[SY_IND] * ds_mat[ind] + geof[TY_IND] * dt_mat[ind]);
        D[ind] += *nz * (geof[RZ_IND] * dr_mat[ind] + geof[SZ_IND] * ds_mat[ind] + geof[TZ_IND] * dt_mat[ind]);
      }
    }

    const DG_FP gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 2) * *fscale;

    DG_FP tmp[DG_NP];
    for(int i = 0; i < DG_NP; i++) {
      tmp[i] = 0.0;
      for(int j = 0; j < DG_NPF; j++) {
        // int mm_ind = i + fmaskB[j] * DG_NP;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], DG_NP, DG_NP);
        u_rhs[i] += gtau * *sJ * mmF[mm_ind] * u_bc[j];
        tmp[i] += *sJ * mmF[mm_ind] * u_bc[j];
      }
    }

    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int mm_ind = i + fmaskB[j] * DG_NP;
        int d_ind = DG_MAT_IND(j, i, DG_NP, DG_NP);
        u_rhs[i] += -D[d_ind] * tmp[j];
      }
    }
  } else if(*u_bc_type == 1) {
    // Neumann
    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NPF; j++) {
        // int mm_ind = i + fmaskB[j] * DG_NP;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], DG_NP, DG_NP);
        u_rhs[i] += *sJ * mmF[mm_ind] * u_bc[j];
      }
    }
  } else {
    // Do nothing for BC_SLIP 
    // not implementing friction, so tangent Neumann portion = 0
    // not porous, so normal Dirichlet portion = 0
  }

  // v bcs
  if(*v_bc_type == 0) {
    // Dirichlet
    DG_FP D[DG_NP * DG_NP];
    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int ind = i + j * DG_NP;
        int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);

        D[ind] = *nx * (geof[RX_IND] * dr_mat[ind] + geof[SX_IND] * ds_mat[ind] + geof[TX_IND] * dt_mat[ind]);
        D[ind] += *ny * (geof[RY_IND] * dr_mat[ind] + geof[SY_IND] * ds_mat[ind] + geof[TY_IND] * dt_mat[ind]);
        D[ind] += *nz * (geof[RZ_IND] * dr_mat[ind] + geof[SZ_IND] * ds_mat[ind] + geof[TZ_IND] * dt_mat[ind]);
      }
    }

    const DG_FP gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 2) * *fscale;

    DG_FP tmp[DG_NP];
    for(int i = 0; i < DG_NP; i++) {
      tmp[i] = 0.0;
      for(int j = 0; j < DG_NPF; j++) {
        // int mm_ind = i + fmaskB[j] * DG_NP;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], DG_NP, DG_NP);
        v_rhs[i] += gtau * *sJ * mmF[mm_ind] * v_bc[j];
        tmp[i] += *sJ * mmF[mm_ind] * v_bc[j];
      }
    }

    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int mm_ind = i + fmaskB[j] * DG_NP;
        int d_ind = DG_MAT_IND(j, i, DG_NP, DG_NP);
        v_rhs[i] += -D[d_ind] * tmp[j];
      }
    }
  } else if(*v_bc_type == 1) {
    // Neumann
    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NPF; j++) {
        // int mm_ind = i + fmaskB[j] * DG_NP;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], DG_NP, DG_NP);
        v_rhs[i] += *sJ * mmF[mm_ind] * v_bc[j];
      }
    }
  } else {
    // Do nothing for BC_SLIP 
    // not implementing friction, so tangent Neumann portion = 0
    // not porous, so normal Dirichlet portion = 0
  }

  // w bcs
  if(*w_bc_type == 0) {
    // Dirichlet
    DG_FP D[DG_NP * DG_NP];
    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int ind = i + j * DG_NP;
        int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);

        D[ind] = *nx * (geof[RX_IND] * dr_mat[ind] + geof[SX_IND] * ds_mat[ind] + geof[TX_IND] * dt_mat[ind]);
        D[ind] += *ny * (geof[RY_IND] * dr_mat[ind] + geof[SY_IND] * ds_mat[ind] + geof[TY_IND] * dt_mat[ind]);
        D[ind] += *nz * (geof[RZ_IND] * dr_mat[ind] + geof[SZ_IND] * ds_mat[ind] + geof[TZ_IND] * dt_mat[ind]);
      }
    }

    const DG_FP gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 2) * *fscale;

    DG_FP tmp[DG_NP];
    for(int i = 0; i < DG_NP; i++) {
      tmp[i] = 0.0;
      for(int j = 0; j < DG_NPF; j++) {
        // int mm_ind = i + fmaskB[j] * DG_NP;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], DG_NP, DG_NP);
        w_rhs[i] += gtau * *sJ * mmF[mm_ind] * w_bc[j];
        tmp[i] += *sJ * mmF[mm_ind] * w_bc[j];
      }
    }

    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int mm_ind = i + fmaskB[j] * DG_NP;
        int d_ind = DG_MAT_IND(j, i, DG_NP, DG_NP);
        w_rhs[i] += -D[d_ind] * tmp[j];
      }
    }
  } else if(*w_bc_type == 1) {
    // Neumann
    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NPF; j++) {
        // int mm_ind = i + fmaskB[j] * DG_NP;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], DG_NP, DG_NP);
        w_rhs[i] += *sJ * mmF[mm_ind] * w_bc[j];
      }
    }
  } else {
    // Do nothing for BC_SLIP 
    // not implementing friction, so tangent Neumann portion = 0
    // not porous, so normal Dirichlet portion = 0
  }
}