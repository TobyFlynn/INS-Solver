inline void ins_pressure_bc_2d(const DG_FP *t, const int *bedge_type,
                               const int *faceNum, const DG_FP *nx,
                               const DG_FP *ny, const DG_FP *fscale, const DG_FP *x,
                               const DG_FP *y, const DG_FP *N0,
                               const DG_FP *N1, const DG_FP *gradCurlVel0,
                               const DG_FP *gradCurlVel1, DG_FP *dPdN) {
  // Get constants for this element's order
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];

  if(*bedge_type == BC_TYPE_NO_SLIP || *bedge_type == BC_TYPE_SLIP_X || *bedge_type == BC_TYPE_SLIP_Y) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      const int find = *faceNum * DG_NPF + i;
      DG_FP res1 = -N0[fmask_ind] - gradCurlVel1[fmask_ind] / r_ynolds;
      DG_FP res2 = -N1[fmask_ind] + gradCurlVel0[fmask_ind] / r_ynolds;
      dPdN[find] += *fscale * (*nx * res1 + *ny * res2);
    }
  } else {
    int tmp_bc_type = -1;
    ps2d_custom_bc_get_pr_type(*bedge_type, tmp_bc_type);
    if(tmp_bc_type == BC_NEUMANN) {
      for(int i = 0; i < DG_NPF; i++) {
        const int fmask_ind = fmask[i];
        DG_FP bc_neumann = ps2d_custom_bc_get_pr_neumann(*bedge_type, *t, x[fmask_ind], y[fmask_ind], *nx, *ny,
                                N0[fmask_ind], N1[fmask_ind], gradCurlVel0[fmask_ind], gradCurlVel1[fmask_ind], r_ynolds);
        const int find = *faceNum * DG_NPF + i;
        dPdN[find] += *fscale * bc_neumann;
      }
    }
  }
}
