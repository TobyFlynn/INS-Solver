inline void vmf_2d_apply_bc(const int *p, const int *faceNum, const int *u_bc_type,
                            const int *v_bc_type, const DG_FP *nx, const DG_FP *ny, 
                            const DG_FP *fscale, const DG_FP *sJ, const DG_FP *geof, 
                            const DG_FP *u_bc, const DG_FP *v_bc, DG_FP *u_rhs, 
                            DG_FP *v_rhs) {
  const DG_FP *dr_mat = &dg_Dr_kernel[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &dg_Ds_kernel[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *mmF0_mat = &dg_MM_F0_kernel[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *mmF1_mat = &dg_MM_F1_kernel[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *mmF2_mat = &dg_MM_F2_kernel[(*p - 1) * DG_NP * DG_NP];
  const int dg_np  = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 1];

  const DG_FP *mmF;
  if(*faceNum == 0)
    mmF = mmF0_mat;
  else if(*faceNum == 1)
    mmF = mmF1_mat;
  else
    mmF = mmF2_mat;

  const int find = *faceNum * dg_npf;
  const int *fmask  = &FMASK[(*p - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * dg_npf];

  if(*u_bc_type == 0) {
    // Dirichlet
    const DG_FP rx = geof[RX_IND];
    const DG_FP ry = geof[RY_IND];
    const DG_FP sx = geof[SX_IND];
    const DG_FP sy = geof[SY_IND];
    DG_FP D[DG_NP * DG_NP];
    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_np; j++) {
        // int ind = i + j * dg_np;
        int ind = DG_MAT_IND(i, j, dg_np, dg_np);

        D[ind] = *nx * (rx * dr_mat[ind] + sx * ds_mat[ind]);
        D[ind] += *ny * (ry * dr_mat[ind] + sy * ds_mat[ind]);
      }
    }

    const DG_FP gtau = 2.0 * (*p + 1) * (*p + 2) * *fscale;

    DG_FP tmp[DG_NP];
    for(int i = 0; i < dg_np; i++) {
      tmp[i] = 0.0;
      for(int j = 0; j < dg_npf; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], dg_np, dg_np);
        u_rhs[i] += gtau * *sJ * mmF[mm_ind] * u_bc[j];
        tmp[i] += *sJ * mmF[mm_ind] * u_bc[j];
      }
    }

    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_np; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int d_ind = DG_MAT_IND(j, i, dg_np, dg_np);
        u_rhs[i] += -D[d_ind] * tmp[j];
      }
    }
  } else if(*u_bc_type == 1) {
    // Neumann
    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_npf; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], dg_np, dg_np);
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
    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_np; j++) {
        // int ind = i + j * dg_np;
        int ind = DG_MAT_IND(i, j, dg_np, dg_np);

        D[ind] = *nx * (rx * dr_mat[ind] + sx * ds_mat[ind]);
        D[ind] += *ny * (ry * dr_mat[ind] + sy * ds_mat[ind]);
      }
    }

    const DG_FP gtau = 2.0 * (*p + 1) * (*p + 2) * *fscale;

    DG_FP tmp[DG_NP];
    for(int i = 0; i < dg_np; i++) {
      tmp[i] = 0.0;
      for(int j = 0; j < dg_npf; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], dg_np, dg_np);
        v_rhs[i] += gtau * *sJ * mmF[mm_ind] * v_bc[j];
        tmp[i] += *sJ * mmF[mm_ind] * v_bc[j];
      }
    }

    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_np; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int d_ind = DG_MAT_IND(j, i, dg_np, dg_np);
        v_rhs[i] += -D[d_ind] * tmp[j];
      }
    }
  } else if(*v_bc_type == 1) {
    // Neumann
    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_npf; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], dg_np, dg_np);
        v_rhs[i] += *sJ * mmF[mm_ind] * v_bc[j];
      }
    }
  } else {
    // Do nothing for BC_SLIP 
    // not implementing friction, so tangent Neumann portion = 0
    // not porous, so normal Dirichlet portion = 0
  }
}
