inline void temperature_diff_2d_bc_0(const int *bc_type, const int *edgeNum,
                      const DG_FP *nx, const DG_FP *ny, const DG_FP *fscale,
                      const DG_FP *x, const DG_FP *y, const DG_FP *val,
                      DG_FP *flux_x, DG_FP *flux_y) {
  // Work out which edge for each element
  const int edge = *edgeNum;
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];

  DG_FP vP[DG_NPF];
  if(*bc_type == BC_TYPE_NO_SLIP || *bc_type == BC_TYPE_SLIP || *bc_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      vP[i] = val[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      vP[i] = ps2d_custom_bc_get_temperature(*bc_type, x[fmask_ind], y[fmask_ind], val[fmask_ind]);
    }
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    const int out_ind = edge * DG_NPF + i;

    flux_x[out_ind] = *fscale * *nx * 0.5 * (val[fmask_ind] + vP[i]);
    flux_y[out_ind] = *fscale * *ny * 0.5 * (val[fmask_ind] + vP[i]);
  }
}
