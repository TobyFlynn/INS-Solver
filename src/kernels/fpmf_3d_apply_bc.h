inline void fpmf_3d_apply_bc(const int *p, const DG_FP *dr, const DG_FP *ds, const DG_FP *dt,
                             const DG_FP *mmF0, const DG_FP *mmF1, const DG_FP *mmF2,
                             const DG_FP *mmF3, const int *faceNum, const int *bc_type,
                             const DG_FP *nx, const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                             const DG_FP *sJ, const DG_FP *rx, const DG_FP *sx, const DG_FP *tx,
                             const DG_FP *ry, const DG_FP *sy, const DG_FP *ty, const DG_FP *rz,
                             const DG_FP *sz, const DG_FP *tz, const DG_FP *fact, const DG_FP *bc,
                             DG_FP *rhs) {
  const DG_FP *dr_mat = &dr[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dt[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *mmF0_mat = &mmF0[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *mmF1_mat = &mmF1[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *mmF2_mat = &mmF2[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *mmF3_mat = &mmF3[(*p - 1) * DG_NP * DG_NP];
  const int dg_np  = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 1];

  const DG_FP tau_order = (DG_FP) *p; // (DG_FP) DG_ORDER;

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
  const int *fmask  = &FMASK[(*p - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * dg_npf];

  if(*bc_type == 0) {
    // Dirichlet
    DG_FP D[DG_NP * DG_NP];
    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_np; j++) {
        // int ind = i + j * dg_np;
        int ind = DG_MAT_IND(i, j, dg_np, dg_np);

        D[ind] = *nx * (rx[0] * dr_mat[ind] + sx[0] * ds_mat[ind] + tx[0] * dt_mat[ind]);
        D[ind] += *ny * (ry[0] * dr_mat[ind] + sy[0] * ds_mat[ind] + ty[0] * dt_mat[ind]);
        D[ind] += *nz * (rz[0] * dr_mat[ind] + sz[0] * ds_mat[ind] + tz[0] * dt_mat[ind]);
        D[ind] *= fact[i];
      }
    }

    DG_FP gtau = fact[fmaskB[0]];
    for(int i = 1; i < dg_npf; i++) {
      gtau = fmax(gtau, fact[fmaskB[i]]);
    }
    gtau *= 2.0 * (tau_order + 1) * (tau_order + 2) * *fscale;

    DG_FP tmp[DG_NP];
    for(int i = 0; i < dg_np; i++) {
      tmp[i] = 0.0;
      for(int j = 0; j < dg_npf; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], dg_np, dg_np);
        rhs[i] += gtau * *sJ * mmF[mm_ind] * bc[j];
        tmp[i] += *sJ * mmF[mm_ind] * bc[j];
      }
    }

    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_np; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int d_ind = DG_MAT_IND(j, i, dg_np, dg_np);
        rhs[i] += -D[d_ind] * tmp[j];
      }
    }
  } else {
    // Neumann
    for(int i = 0; i < dg_np; i++) {
      for(int j = 0; j < dg_npf; j++) {
        // int mm_ind = i + fmaskB[j] * dg_np;
        int mm_ind = DG_MAT_IND(i, fmaskB[j], dg_np, dg_np);
        rhs[i] += *sJ * mmF[mm_ind] * bc[j];
      }
    }
  }
}
