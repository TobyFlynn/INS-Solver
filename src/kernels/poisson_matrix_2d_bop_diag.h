inline void poisson_matrix_2d_bop_diag(const int *order, const int *faceNum,
                                       const int *bc_type, const DG_FP *nx,
                                       const DG_FP *ny, const DG_FP *fscale,
                                       const DG_FP *sJ, const DG_FP *geof,
                                       DG_FP *diag) {
  // Do nothing for Neumann boundary conditions
  if(*bc_type == 1)
    return;

  // Handle Dirichlet boundary conditions
  const DG_FP *dr_mat = &dg_Dr_kernel[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &dg_Ds_kernel[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *mmF0_mat = &dg_MM_F0_kernel[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *mmF1_mat = &dg_MM_F1_kernel[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *mmF2_mat = &dg_MM_F2_kernel[(*order - 1) * DG_NP * DG_NP];
  const int dg_np  = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS + 1];

  const DG_FP *mmF;
  if(*faceNum == 0)
    mmF = mmF0_mat;
  else if(*faceNum == 1)
    mmF = mmF1_mat;
  else
    mmF = mmF2_mat;

  const int find = *faceNum * dg_npf;
  const int *fmask  = &FMASK[(*order - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * dg_npf];

  DG_FP D[DG_NP * DG_NP];
  const DG_FP r_fact = nx[0] * geof[RX_IND] + ny[0] * geof[RY_IND];
  const DG_FP s_fact = nx[0] * geof[SX_IND] + ny[0] * geof[SY_IND];
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      D[ind] = r_fact * dr_mat[ind] + s_fact * ds_mat[ind];
    }
  }

  const DG_FP gtau = 2.0 * (*order + 1) * (*order + 2) * *fscale;

  for(int i = 0; i < dg_np; i++) {
    DG_FP tmp = 0.0;
    for(int k = 0; k < dg_np; k++) {
      // int a_ind0 = i + k * dg_np;
      int a_ind0 = DG_MAT_IND(i, k, dg_np, dg_np);
      // int a_ind1 = i * dg_np + k;
      int a_ind1 = DG_MAT_IND(k, i, dg_np, dg_np);
      // int b_ind  = j * dg_np + k;
      int b_ind = DG_MAT_IND(k, i, dg_np, dg_np);
      tmp += mmF[a_ind0] * D[b_ind] + D[a_ind1] * mmF[b_ind];
    }
    int op_ind = DG_MAT_IND(i, i, dg_np, dg_np);
    diag[i] += gtau * *sJ * mmF[op_ind] - *sJ * tmp;
  }
}
