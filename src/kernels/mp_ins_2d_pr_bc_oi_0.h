inline void mp_ins_2d_pr_bc_oi_0(const DG_FP *t, const int *bedge_type,
                           const int *faceNum, const DG_FP *nx, const DG_FP *ny,
                           const DG_FP *fscale, const DG_FP *x, const DG_FP *y,
                           const DG_FP *N0, const DG_FP *N1,
                           const DG_FP *gradCurlVel0, const DG_FP *gradCurlVel1,
                           const DG_FP *rho, DG_FP *dPdN) {
  const int edge = *faceNum;
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];
  const int fInd = edge * DG_NPF;
  const int fIndCub = edge * DG_CUB_SURF_2D_NP;

  DG_FP oi_N0[DG_CUB_SURF_2D_NP], oi_N1[DG_CUB_SURF_2D_NP];
  DG_FP oi_gradCurlVel0[DG_CUB_SURF_2D_NP], oi_gradCurlVel1[DG_CUB_SURF_2D_NP];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    oi_N0[i] = 0.0; oi_N1[i] = 0.0;
    oi_gradCurlVel0[i] = 0.0; oi_gradCurlVel1[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];

    const DG_FP _N0 = N0[fmask_ind];
    const DG_FP _N1 = N1[fmask_ind];
    const DG_FP _gradCurlVel0 = gradCurlVel0[fmask_ind];
    const DG_FP _gradCurlVel1 = gradCurlVel1[fmask_ind];

    for(int j = 0; j < DG_CUB_SURF_2D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCub + j, fInd + i, DG_CUB_SURF_2D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf2d_Interp_kernel[ind];
      oi_N0[j] += mat_val * _N0;
      oi_N1[j] += mat_val * _N1;
      oi_gradCurlVel0[j] += mat_val * _gradCurlVel0;
      oi_gradCurlVel1[j] += mat_val * _gradCurlVel1;
    }
  }

  if(*bedge_type == BC_TYPE_NO_SLIP || *bedge_type == BC_TYPE_SLIP || *bedge_type == BC_TYPE_SLIP_X || *bedge_type == BC_TYPE_SLIP_Y) {
    for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
      DG_FP res1 = -oi_N0[i] - oi_gradCurlVel1[i] / (r_ynolds * rho[fIndCub + i]);
      DG_FP res2 = -oi_N1[i] + oi_gradCurlVel0[i] / (r_ynolds * rho[fIndCub + i]);
      dPdN[fIndCub + i] += *fscale * (*nx * res1 + *ny * res2);
    }
  } else {
    int tmp_bc_type = -1;
    ps2d_custom_bc_get_pr_type(*bedge_type, tmp_bc_type);
    if(tmp_bc_type == BC_NEUMANN) {
      // Also need x and y at over int points
      DG_FP oi_x[DG_CUB_SURF_2D_NP], oi_y[DG_CUB_SURF_2D_NP];
      for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
        oi_x[i] = 0.0; oi_y[i] = 0.0;
      }
      for(int i = 0; i < DG_NPF; i++) {
        const int fmask_ind = fmask[i];
        const DG_FP _x = x[fmask_ind];
        const DG_FP _y = y[fmask_ind];
        for(int j = 0; j < DG_CUB_SURF_2D_NP; j++) {
          const int ind = DG_MAT_IND(fIndCub + j, fInd + i, DG_CUB_SURF_2D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
          const DG_FP mat_val = dg_cubSurf2d_Interp_kernel[ind];
          oi_x[j] += mat_val * _x;
          oi_y[j] += mat_val * _y;
        }
      }

      for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
        DG_FP bc_neumann = ps2d_custom_bc_get_pr_neumann_multiphase(*bedge_type, *t, oi_x[i], oi_y[i], *nx, *ny,
                              oi_N0[i], oi_N1[i], oi_gradCurlVel0[i], oi_gradCurlVel1[i], r_ynolds, rho[fIndCub + i]);
        dPdN[fIndCub + i] += *fscale * bc_neumann;
      }
    }
  }
}
