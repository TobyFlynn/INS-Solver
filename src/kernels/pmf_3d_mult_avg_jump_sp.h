inline void pmf_3d_mult_avg_jump_sp(const int *order, const int *bc_type,
                               const int *faceNum, const float *in,
                               const float *in_x, const float *in_y,
                               const float *in_z, float *jump, float *avg_x,
                               float *avg_y, float *avg_z) {
  if(*bc_type == 1)
    return;

  const int p = order[0];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

  for(int j = 0; j < dg_npf; j++) {
    const int fmaskInd = FMASK[(p - 1) * 4 * DG_NPF + faceNum[0] * dg_npf + j];
    const int ind = faceNum[0] * dg_npf + j;
    jump[ind] += 2.0 * in[fmaskInd];
    avg_x[ind] += in_x[fmaskInd];
    avg_y[ind] += in_y[fmaskInd];
    avg_z[ind] += in_z[fmaskInd];
  }
}
