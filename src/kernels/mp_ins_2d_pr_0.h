inline void mp_ins_2d_pr_0(const DG_FP *t, const int *bedge_type,
                           const int *faceNum, const DG_FP *nx, const DG_FP *ny,
                           const DG_FP *fscale, const DG_FP *x, const DG_FP *y,
                           const DG_FP *N0, const DG_FP *N1,
                           const DG_FP *gradCurlVel0, const DG_FP *gradCurlVel1,
                           const DG_FP *rho, DG_FP *dPdN) {
  // Get constants for this element's order
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  const DG_FP PI = 3.141592653589793238463;

  if(*bedge_type == 0 || *bedge_type == 2) {
    // Inflow or Wall
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      const int find = *faceNum * DG_NPF + i;
      DG_FP res1 = -N0[fmask_ind] - gradCurlVel1[fmask_ind] / (r_ynolds * rho[fmask_ind]);
      DG_FP res2 = -N1[fmask_ind] + gradCurlVel0[fmask_ind] / (r_ynolds * rho[fmask_ind]);
      dPdN[find] += *fscale * (*nx * res1 + *ny * res2);
    }
  }
}
