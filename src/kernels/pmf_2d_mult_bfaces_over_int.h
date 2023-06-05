inline void pmf_2d_mult_bfaces_over_int(const int *p, const int *bc_type,
                  const int *edgeNum, const DG_FP *sJ, const DG_FP *nx,
                  const DG_FP *ny, const DG_FP *h, const DG_FP *in,
                  const DG_FP *in_x, const DG_FP *in_y, DG_FP *out,
                  DG_FP *l_x, DG_FP *l_y) {
  if(*bc_type == 1)
    return;

  // Get constants
  // Using same Gauss points so should be able to replace dg_gf_npL and
  // dg_gf_npR with DG_GF_NP
  const int dg_np      = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf     = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 1];
  const int dg_gf_np   = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 4];
  const DG_FP *gaussW = &gaussW_g[(*p - 1) * DG_GF_NP];
  const int exInd = *edgeNum * DG_GF_NP;

  DG_FP max_hinv = h[*edgeNum * dg_npf];
  DG_FP tau = 0.5 * (DG_ORDER + 1) * (DG_ORDER + 2) * max_hinv;

  for(int i = 0; i < DG_GF_NP; i++) {
    int ind = exInd + i;

    const DG_FP diff_u = in[ind];
    const DG_FP diff_u_x = nx[ind] * in_x[ind];
    const DG_FP diff_u_y = ny[ind] * in_y[ind];
    const DG_FP diff_u_grad = diff_u_x + diff_u_y;

    out[ind] += 0.5 * gaussW[i] * sJ[ind] * (tau * diff_u - diff_u_grad);
    const DG_FP l_tmp = 0.5 * gaussW[i] * sJ[ind] * -diff_u;
    l_x[ind] += nx[ind] * l_tmp;
    l_y[ind] += ny[ind] * l_tmp;
  }
}
