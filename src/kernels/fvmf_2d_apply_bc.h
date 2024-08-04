inline void fvmf_2d_apply_bc(const int *faceNum, const int *u_bc_type,
                             const int *v_bc_type, const DG_FP *nx, const DG_FP *ny, 
                             const DG_FP *fscale, const DG_FP *sJ, const DG_FP *geof, 
                             const DG_FP *factor, const DG_FP *u_bc, const DG_FP *v_bc, 
                             DG_FP *u_rhs, DG_FP *v_rhs) {
  const DG_FP *dr_mat = &dg_Dr_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &dg_Ds_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF0_mat = &dg_MM_F0_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF1_mat = &dg_MM_F1_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *mmF2_mat = &dg_MM_F2_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];

  const DG_FP *mmF;
  if(*faceNum == 0)
    mmF = mmF0_mat;
  else if(*faceNum == 1)
    mmF = mmF1_mat;
  else
    mmF = mmF2_mat;

  const int find = *faceNum * DG_NPF;
  const int *fmask  = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * DG_NPF];

  if(*u_bc_type == 0) {
    // Dirichlet
    const DG_FP rx = geof[RX_IND];
    const DG_FP ry = geof[RY_IND];
    const DG_FP sx = geof[SX_IND];
    const DG_FP sy = geof[SY_IND];
    DG_FP D[DG_NP * DG_NP];
    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int ind = i + j * dg_np;
        int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);

        D[ind] = *nx * (rx * dr_mat[ind] + sx * ds_mat[ind]);
        D[ind] += *ny * (ry * dr_mat[ind] + sy * ds_mat[ind]);
        D[ind] *= factor[i];
      }
    }

    const int fmask_ind_0 = fmaskB[0];
    DG_FP gtau = factor[fmask_ind_0];
    for(int i = 1; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      gtau = fmax(gtau, factor[fmask_ind]);
    }
    gtau *= 2.0 * (DG_ORDER + 1) * (DG_ORDER + 2) * *fscale;

    DG_FP tmp[DG_NP];
    for(int i = 0; i < DG_NP; i++) {
      tmp[i] = 0.0;
      for(int j = 0; j < DG_NPF; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], DG_NP, DG_NP);
        u_rhs[i] += gtau * *sJ * mmF[mm_ind] * u_bc[j];
        tmp[i] += *sJ * mmF[mm_ind] * u_bc[j];
      }
    }

    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int d_ind = DG_MAT_IND(j, i, DG_NP, DG_NP);
        u_rhs[i] += -D[d_ind] * tmp[j];
      }
    }
  } else if(*u_bc_type == 1) {
    // Neumann
    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NPF; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], DG_NP, DG_NP);
        u_rhs[i] += *sJ * mmF[mm_ind] * u_bc[j];
      }
    }
  } else {
    // Do nothing for BC_SLIP 
    // not implementing friction, so tangent Neumann portion = 0
    // not porous, so normal Dirichlet portion = 0
  }

  if(*v_bc_type == 0) {
    // Dirichlet
    const DG_FP rx = geof[RX_IND];
    const DG_FP ry = geof[RY_IND];
    const DG_FP sx = geof[SX_IND];
    const DG_FP sy = geof[SY_IND];
    DG_FP D[DG_NP * DG_NP];
    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int ind = i + j * dg_np;
        int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);

        D[ind] = *nx * (rx * dr_mat[ind] + sx * ds_mat[ind]);
        D[ind] += *ny * (ry * dr_mat[ind] + sy * ds_mat[ind]);
        D[ind] *= factor[i];
      }
    }

    const int fmask_ind_0 = fmaskB[0];
    DG_FP gtau = factor[fmask_ind_0];
    for(int i = 1; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      gtau = fmax(gtau, factor[fmask_ind]);
    }
    gtau *= 2.0 * (DG_ORDER + 1) * (DG_ORDER + 2) * *fscale;

    DG_FP tmp[DG_NP];
    for(int i = 0; i < DG_NP; i++) {
      tmp[i] = 0.0;
      for(int j = 0; j < DG_NPF; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], DG_NP, DG_NP);
        v_rhs[i] += gtau * *sJ * mmF[mm_ind] * v_bc[j];
        tmp[i] += *sJ * mmF[mm_ind] * v_bc[j];
      }
    }

    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NP; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int d_ind = DG_MAT_IND(j, i, DG_NP, DG_NP);
        v_rhs[i] += -D[d_ind] * tmp[j];
      }
    }
  } else if(*v_bc_type == 1) {
    // Neumann
    for(int i = 0; i < DG_NP; i++) {
      for(int j = 0; j < DG_NPF; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], DG_NP, DG_NP);
        v_rhs[i] += *sJ * mmF[mm_ind] * v_bc[j];
      }
    }
  } else {
    // Do nothing for BC_SLIP 
    // not implementing friction, so tangent Neumann portion = 0
    // not porous, so normal Dirichlet portion = 0
  }
}
