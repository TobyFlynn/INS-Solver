inline void temperature_diff_2d_bc_1(const int *bc_type, const int *edgeNum,
                      const DG_FP *nx, const DG_FP *ny, const DG_FP *fscale,
                      const DG_FP *x, const DG_FP *y, const DG_FP *val_x,
                      const DG_FP *val_y, const DG_FP *val, const DG_FP *vis,
                      DG_FP *flux) {
  // Work out which edge for each element
  const int edge = *edgeNum;
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];

  DG_FP vP[DG_NPF], vPx[DG_NPF], vPy[DG_NPF];
  if(*bc_type == BC_TYPE_NO_SLIP || *bc_type == BC_TYPE_SLIP || *bc_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      vP[i] = val[fmask_ind];
      vPx[i] = val_x[fmask_ind];
      vPy[i] = val_y[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      vP[i] = ps2d_custom_bc_get_temperature(*bc_type, x[fmask_ind], y[fmask_ind], val[fmask_ind]);
      ps2d_custom_bc_get_temperature_grad(*bc_type, x[fmask_ind], y[fmask_ind], val_x[fmask_ind], val_y[fmask_ind], vPx[i], vPy[i]);
    }
  }

  const DG_FP pen = DG_ORDER * DG_ORDER * *fscale;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    const int out_ind = edge * DG_NPF + i;

    const DG_FP avg_val_x = 0.5 * (val_x[fmask_ind] + vPx[i]);
    const DG_FP avg_val_y = 0.5 * (val_y[fmask_ind] + vPy[i]);
    const DG_FP diff_val_0 = val[fmask_ind] - vP[i];

    flux[out_ind] = *fscale * (*nx * avg_val_x + *ny * avg_val_y - pen * vis[fmask_ind] * diff_val_0);
  }
}
