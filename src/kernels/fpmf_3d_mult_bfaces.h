inline void fpmf_3d_mult_bfaces(const int *order, const int *bc_type, const int *faceNum,
                                const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                                const DG_FP *fscale, const DG_FP *sJ, const DG_FP *in,
                                const DG_FP *fact, const DG_FP *in_x, const DG_FP *in_y,
                                const DG_FP *in_z, DG_FP *l_x, DG_FP *l_y, DG_FP *l_z,
                                DG_FP *out) {
  if(*bc_type == 1)
    return;

  const int p = order[0];
  const int dg_np  = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

  const DG_FP tau_order = (DG_FP) p; // (DG_FP) DG_ORDER;

  const int find = faceNum[0] * dg_npf;
  const int *fmask  = &FMASK[(p - 1) * 4 * DG_NPF];
  const int *fmaskB = &fmask[faceNum[0] * dg_npf];

  const int fmask_ind_0 = fmaskB[0];
  DG_FP gtau = fact[fmask_ind_0];
  for(int i = 1; i < dg_npf; i++) {
    const int fmask_ind = fmaskB[i];
    gtau = fmax(gtau, fact[fmask_ind]);
  }
  gtau *= 2.0 * (tau_order + 1) * (tau_order + 2) * *fscale;

  for(int j = 0; j < dg_npf; j++) {
    const int fmaskInd = fmaskB[j];
    const DG_FP diff_u = in[fmaskInd];
    const DG_FP diff_u_x = nx[0] * in_x[fmaskInd];
    const DG_FP diff_u_y = ny[0] * in_y[fmaskInd];
    const DG_FP diff_u_z = nz[0] * in_z[fmaskInd];
    const DG_FP diff_u_grad = diff_u_x + diff_u_y + diff_u_z;

    const int ind = find + j;
    out[ind] += sJ[0] * (gtau * diff_u - diff_u_grad);
    const DG_FP l_tmp = fact[fmaskInd] * sJ[0] * -diff_u;
    l_x[ind] += nx[0] * l_tmp;
    l_y[ind] += ny[0] * l_tmp;
    l_z[ind] += nz[0] * l_tmp;
  }
}
