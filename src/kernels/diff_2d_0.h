inline void diff_2d_0(const int *edgeNum, const bool *rev, const DG_FP *nx,
                      const DG_FP *ny, const DG_FP *fscale, const DG_FP **val, 
                      DG_FP **flux_x, DG_FP **flux_y) {
  // Work out which edge for each element
  const int edgeL = edgeNum[0];
  const int edgeR = edgeNum[1];
  const bool reverse = *rev;
  const int *fmaskL = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeL * DG_NPF];
  const int *fmaskR = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeR * DG_NPF];

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = reverse ? fmaskR[DG_NPF - i - 1] : fmaskR[i];
    const int outL_ind = edgeL * DG_NPF + i;
    const int outR_ind = reverse ? edgeR * DG_NPF + DG_NPF - i - 1 : edgeR * DG_NPF + i;

    DG_FP avg_val = 0.5 * (val[0][fmaskL_ind] + val[1][fmaskR_ind]);

    flux_x[0][outL_ind] = fscale[0] * nx[0] * avg_val;
    flux_x[1][outR_ind] = fscale[1] * nx[1] * avg_val;
    flux_y[0][outL_ind] = fscale[0] * ny[0] * avg_val;
    flux_y[1][outR_ind] = fscale[1] * ny[1] * avg_val;
  }
}
