inline void poisson_coarse_matrix_3d_bop(const DG_FP *dr,
                                         const DG_FP *ds, const DG_FP *dt,
                                         const DG_FP *mmF0, const DG_FP *mmF1,
                                         const DG_FP *mmF2, const DG_FP *mmF3,
                                         const int *faceNum, const int *bc_type,
                                         const DG_FP *nx, const DG_FP *ny,
                                         const DG_FP *nz, const DG_FP *fscale,
                                         const DG_FP *sJ, const DG_FP *rx,
                                         const DG_FP *sx, const DG_FP *tx,
                                         const DG_FP *ry, const DG_FP *sy,
                                         const DG_FP *ty, const DG_FP *rz,
                                         const DG_FP *sz, const DG_FP *tz,
                                         DG_FP *op1) {
  // Do nothing for Neumann boundary conditions
  if(*bc_type == 1)
    return;

  // Handle Dirichlet boundary conditions
  const DG_FP *dr_mat = dr;
  const DG_FP *ds_mat = ds;
  const DG_FP *dt_mat = dt;
  const DG_FP *mmF0_mat = mmF0;
  const DG_FP *mmF1_mat = mmF1;
  const DG_FP *mmF2_mat = mmF2;
  const DG_FP *mmF3_mat = mmF3;

  const DG_FP *mmF;
  if(*faceNum == 0)
    mmF = mmF0_mat;
  else if(*faceNum == 1)
    mmF = mmF1_mat;
  else if(*faceNum == 2)
    mmF = mmF2_mat;
  else
    mmF = mmF3_mat;

  const int find = *faceNum * DG_NPF_N1;
  const int *fmask  = FMASK;
  const int *fmaskB = &fmask[*faceNum * DG_NPF_N1];

  DG_FP D[DG_NP_N1 * DG_NP_N1];
  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int ind = i + j * DG_NP_N1;
      int ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);

      D[ind] = *nx * (rx[0] * dr_mat[ind] + sx[0] * ds_mat[ind] + tx[0] * dt_mat[ind]);
      D[ind] += *ny * (ry[0] * dr_mat[ind] + sy[0] * ds_mat[ind] + ty[0] * dt_mat[ind]);
      D[ind] += *nz * (rz[0] * dr_mat[ind] + sz[0] * ds_mat[ind] + tz[0] * dt_mat[ind]);
    }
  }

  const DG_FP gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 2) * *fscale;

  for(int i = 0; i < DG_NP_N1; i++) {
    for(int j = 0; j < DG_NP_N1; j++) {
      // int op_ind = i + j * DG_NP_N1;
      int op_ind = DG_MAT_IND(i, j, DG_NP_N1, DG_NP_N1);
      DG_FP tmp = 0.0;
      for(int k = 0; k < DG_NP_N1; k++) {
        // int a_ind0 = i + k * DG_NP_N1;
        int a_ind0 = DG_MAT_IND(i, k, DG_NP_N1, DG_NP_N1);
        // int a_ind1 = i * DG_NP_N1 + k;
        int a_ind1 = DG_MAT_IND(k, i, DG_NP_N1, DG_NP_N1);
        // int b_ind  = j * DG_NP_N1 + k;
        int b_ind = DG_MAT_IND(k, j, DG_NP_N1, DG_NP_N1);
        tmp += -*sJ * mmF[a_ind0] * D[b_ind] - D[a_ind1] * *sJ * mmF[b_ind];
      }
      op1[op_ind] += gtau * *sJ * mmF[op_ind] + tmp;
    }
  }
}
