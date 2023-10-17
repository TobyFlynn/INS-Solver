inline void diff_2d_1(const int *edgeNum, const DG_FP *nx, const DG_FP *ny, 
                      const DG_FP *fscale, const DG_FP *val, DG_FP *flux_x, 
                      DG_FP *flux_y) {
  // Work out which edge for each element
  const int edge = *edgeNum;
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    const int out_ind = edge * DG_NPF + i;

    flux_x[out_ind] = *fscale * *nx * val[fmask_ind];
    flux_y[out_ind] = *fscale * *ny * val[fmask_ind];
  }
}
