inline void ins_3d_pr_2(const DG_FP *t, const int *bc_type, const int *faceNum,
                        const DG_FP *x, const DG_FP *y, const DG_FP *z,
                        int *type, DG_FP *bcs) {
  for(int i = 0; i < DG_NPF; i++) {
    bcs[i] = 0.0;
  }

  if(*bc_type == BC_TYPE_NO_SLIP || *bc_type == BC_TYPE_SLIP) {
    *type = BC_NEUMANN;
  } else if(*bc_type == BC_TYPE_NATURAL_OUTFLOW) {
    *type = BC_DIRICHLET;
  } else {
    ps3d_custom_bc_get_pr_type(*bc_type, *type);

    if(*type == BC_DIRICHLET) {
      const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
      for(int i = 0; i < DG_NPF; i++) {
        const int fmask_ind = fmask[i];
        bcs[i] = ps3d_custom_bc_get_pr_dirichlet(*bc_type, *t, x[fmask_ind], y[fmask_ind], z[fmask_ind]);
      }
    }
  }
}
