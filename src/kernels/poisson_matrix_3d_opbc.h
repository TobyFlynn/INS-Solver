inline void poisson_matrix_3d_opbc(const int *order, const DG_FP *dr,
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
                                   DG_FP *op) {
  const DG_FP *dr_mat = &dr[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dt[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *mmF0_mat = &mmF0[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *mmF1_mat = &mmF1[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *mmF2_mat = &mmF2[(*order - 1) * DG_NP * DG_NP];
  const DG_FP *mmF3_mat = &mmF3[(*order - 1) * DG_NP * DG_NP];
  const int dg_np  = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS + 1];

  const DG_FP *mmF;
  if(*faceNum == 0)
    mmF = mmF0_mat;
  else if(*faceNum == 1)
    mmF = mmF1_mat;
  else if(*faceNum == 2)
    mmF = mmF2_mat;
  else
    mmF = mmF3_mat;

  const int find = *faceNum * dg_npf;
  const int *fmask  = &FMASK[(*order - 1) * 4 * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * dg_npf];

  if(*bc_type == 0) {
    // Dirichlet
    DG_FP D[DG_NP * DG_NP];
    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_np; j++) {
        int ind  = i + j * dg_np;

        D[ind] = *nx * (rx[0] * dr_mat[ind] + sx[0] * ds_mat[ind] + tx[0] * dt_mat[ind]);
        D[ind] += *ny * (ry[0] * dr_mat[ind] + sy[0] * ds_mat[ind] + ty[0] * dt_mat[ind]);
        D[ind] += *nz * (rz[0] * dr_mat[ind] + sz[0] * ds_mat[ind] + tz[0] * dt_mat[ind]);
      }
    }

    const DG_FP gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 1) * *fscale;

    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_npf; j++) {
        int op_ind = i + j * dg_np;
        int mm_ind = i + fmaskB[j] * dg_np;
        op[op_ind] = gtau * *sJ * mmF[mm_ind];
      }
    }

    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_npf; j++) {
        int op_ind = i + j * dg_np;
        for(int k = 0; k < dg_np; k++) {
          int a_ind = i * dg_np + k;
          int b_ind  = fmaskB[j] * dg_np + k;
          op[op_ind] += -D[a_ind] * *sJ * mmF[b_ind];
        }
      }
    }
  } else {
    // Neumann
    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_npf; j++) {
        int op_ind = i + j * dg_np;
        int mm_ind = i + fmaskB[j] * dg_np;
        op[op_ind] = *sJ * mmF[mm_ind];
      }
    }
  }
}
