inline void diff_2d_2(const int *edgeNum, const bool *rev, const DG_FP *nx,
                      const DG_FP *ny, const DG_FP *fscale, const DG_FP **val_x, 
                      const DG_FP **val_y, const DG_FP **val, const DG_FP **vis,
                      DG_FP **flux) {
  // Work out which edge for each element
  const int edgeL = edgeNum[0];
  const int edgeR = edgeNum[1];
  const bool reverse = *rev;
  const int *fmaskL = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeL * DG_NPF];
  const int *fmaskR = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeR * DG_NPF];

  const DG_FP pen = DG_ORDER * DG_ORDER * fmax(fscale[0], fscale[1]);

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = reverse ? fmaskR[DG_NPF - i - 1] : fmaskR[i];
    const int outL_ind = edgeL * DG_NPF + i;
    const int outR_ind = reverse ? edgeR * DG_NPF + DG_NPF - i - 1 : edgeR * DG_NPF + i;

    const DG_FP avg_val_x = 0.5 * (val_x[0][fmaskL_ind] + val_x[1][fmaskR_ind]);
    const DG_FP avg_val_y = 0.5 * (val_y[0][fmaskL_ind] + val_y[1][fmaskR_ind]);
    const DG_FP diff_val_0 = val[0][fmaskL_ind] - val[1][fmaskR_ind];
    const DG_FP diff_val_1 = val[1][fmaskR_ind] - val[0][fmaskL_ind];

    flux[0][outL_ind] = fscale[0] * (nx[0] * avg_val_x + ny[0] * avg_val_y - pen * vis[0][fmaskL_ind] * diff_val_0);
    flux[1][outR_ind] = fscale[1] * (nx[1] * avg_val_x + ny[1] * avg_val_y - pen * vis[1][fmaskR_ind] * diff_val_1);
  }
}
