inline void fpmf_2d_calc_tau_bfaces(const int *order, const int *faceNum,
                                    const DG_FP *fscale, const DG_FP *factor,
                                    DG_FP *tau) {
  const int p = order[0];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];
  const int _faceNum = faceNum[0];
  const int fmask_ind_0 = FMASK[(p - 1) * DG_NUM_FACES * DG_NPF + _faceNum * dg_npf];
  DG_FP gtau = factor[fmask_ind_0];
  for(int i = 1; i < dg_npf; i++) {
    const int fmask_ind = FMASK[(p - 1) * DG_NUM_FACES * DG_NPF + _faceNum * dg_npf + i];
    gtau = fmax(gtau, factor[fmask_ind]);
  }
  gtau *= 2.0 * (p + 1) * (p + 2) * fscale[0];
  tau[_faceNum] = gtau;
}
