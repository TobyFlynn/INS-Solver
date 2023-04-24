inline void ins_3d_pr_0(const DG_FP *t, const int *bc_type, const int *faceNum,
                        const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                        const DG_FP *fscale, const DG_FP *x, const DG_FP *y,
                        const DG_FP *z, const DG_FP *n0, const DG_FP *n1,
                        const DG_FP *n2, const DG_FP *curl20,
                        const DG_FP *curl21, const DG_FP *curl22, DG_FP *dPdN) {
  // TODO handle different boundary conditions
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;
  const DG_FP PI = 3.141592653589793238463;

  if(*bc_type == LW_INFLOW_BC || *bc_type == LW_NO_SLIP_WALL_BC || *bc_type == LW_SLIP_WALL_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      DG_FP res0 = -n0[fmask_ind] - curl20[fmask_ind] / r_ynolds;
      DG_FP res1 = -n1[fmask_ind] - curl21[fmask_ind] / r_ynolds;
      DG_FP res2 = -n2[fmask_ind] - curl22[fmask_ind] / r_ynolds;

      dPdN[fInd + i] += *fscale * (*nx * res0 + *ny * res1 + *nz * res2);
    }

    if(*bc_type == LW_INFLOW_BC) {
      for(int i = 0; i < DG_NPF; i++) {
        // dPdN[fInd + i] += *fscale * PI * cos(PI * (*t));
      }
    }
  }
}
