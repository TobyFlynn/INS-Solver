inline void mp_ins_3d_pr_0(const DG_FP *t, const int *bc_type, const int *faceNum,
                           const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                           const DG_FP *fscale, const DG_FP *x, const DG_FP *y,
                           const DG_FP *z, const DG_FP *n0, const DG_FP *n1,
                           const DG_FP *n2, const DG_FP *curl20,
                           const DG_FP *curl21, const DG_FP *curl22,
                           const DG_FP *rho, DG_FP *dPdN) {
  // TODO handle different boundary conditions
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;

  if(*bc_type == BC_TYPE_NO_SLIP || *bc_type == BC_TYPE_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      DG_FP res0 = -n0[fmask_ind] - curl20[fmask_ind] / (r_ynolds * rho[fmask_ind]);
      DG_FP res1 = -n1[fmask_ind] - curl21[fmask_ind] / (r_ynolds * rho[fmask_ind]);
      DG_FP res2 = -n2[fmask_ind] - curl22[fmask_ind] / (r_ynolds * rho[fmask_ind]);

      dPdN[fInd + i] += *fscale * (*nx * res0 + *ny * res1 + *nz * res2);
    }
  } else {
    int tmp_bc_type = -1;
    ps3d_custom_bc_get_pr_type(*bc_type, tmp_bc_type);
    if(tmp_bc_type == BC_NEUMANN) {
      for(int i = 0; i < DG_NPF; i++) {
        const int fmask_ind = fmask[i];
        DG_FP bc_neumann = ps3d_custom_bc_get_pr_neumann_multiphase(*bc_type, *t, x[fmask_ind], y[fmask_ind],
                              z[fmask_ind], *nx, *ny, *nz, n0[fmask_ind], n1[fmask_ind], n2[fmask_ind],
                              curl20[fmask_ind], curl21[fmask_ind], curl22[fmask_ind], r_ynolds, rho[fmask_ind]);
        const int find = *faceNum * DG_NPF + i;
        dPdN[find] += *fscale * bc_neumann;
      }
    }
  }
}
