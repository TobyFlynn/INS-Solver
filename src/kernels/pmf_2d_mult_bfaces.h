inline void pmf_2d_mult_bfaces(const int *order, const int *bc_type, const int *faceNum,
                               const DG_FP *nx, const DG_FP *ny, const DG_FP *fscale,
                               const DG_FP *sJ, const DG_FP *in, const DG_FP *in_x,
                               const DG_FP *in_y, DG_FP *l_x, DG_FP *l_y, DG_FP *out) {
  if(*bc_type == 1)
    return;

  const int p = *order;
  const int dg_np  = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

  const int find = faceNum[0] * dg_npf;
  const int *fmask  = &FMASK[(p - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskB = &fmask[faceNum[0] * dg_npf];

  const DG_FP gtau = 2.0 * (p + 1) * (p + 2) * fscale[0];

  for(int j = 0; j < dg_npf; j++) {
    const int fmaskInd = fmaskB[j];
    const DG_FP diff_u = in[fmaskInd];
    const DG_FP diff_u_x = nx[0] * in_x[fmaskInd];
    const DG_FP diff_u_y = ny[0] * in_y[fmaskInd];
    const DG_FP diff_u_grad = diff_u_x + diff_u_y;

    const int ind = find + j;
    out[ind] += sJ[0] * (gtau * diff_u - diff_u_grad);
    const DG_FP l_tmp = sJ[0] * -diff_u;
    l_x[ind] += nx[0] * l_tmp;
    l_y[ind] += ny[0] * l_tmp;
  }
}
