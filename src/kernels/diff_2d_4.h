inline void diff_2d_4(const int *edgeNum, const DG_FP *nx, const DG_FP *ny, 
                      const DG_FP *fscale, const DG_FP *val_x, const DG_FP *val_y, 
                      const DG_FP *val, const DG_FP *vis, DG_FP *flux) {
  // Work out which edge for each element
  const int edge = *edgeNum;
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];

  const DG_FP pen = DG_ORDER * DG_ORDER * *fscale;

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    const int out_ind = edge * DG_NPF + i;

    const DG_FP avg_val_x = val_x[fmask_ind];
    const DG_FP avg_val_y = val_y[fmask_ind];
    const DG_FP diff_val_0 = 0.0;

    flux[out_ind] = *fscale * (*nx * avg_val_x + *ny * avg_val_y - pen * vis[fmask_ind] * diff_val_0);
  }
}
