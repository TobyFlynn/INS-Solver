inline void ins_2d_pr_bc(const DG_FP *t, const int *bc_type, const int *type, 
                         const int *bedgeNum, const DG_FP *x, const DG_FP *y, 
                         DG_FP *bcs) {
  if(*type == BC_DIRICHLET && *bc_type != BC_TYPE_NATURAL_OUTFLOW) {
    const int edge = bedgeNum[0];
    const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      bcs[i] = ps2d_custom_bc_get_pr_dirichlet(*bc_type, *t, x[fmask_ind], y[fmask_ind]);
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      bcs[i] = 0.0;
    }
  }
}
