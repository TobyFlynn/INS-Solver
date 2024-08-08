inline void fvmf_2d_calc_tau_bfaces(const int *faceNum, const DG_FP *fscale, const DG_FP *factor, DG_FP *tau) {
  const int _faceNum = faceNum[0];
  const int fmask_ind_0 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + _faceNum * DG_NPF];
  DG_FP gtau = factor[fmask_ind_0];
  for(int i = 1; i < DG_NPF; i++) {
    const int fmask_ind = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + _faceNum * DG_NPF + i];
    gtau = fmax(gtau, factor[fmask_ind]);
  }
  gtau *= 2.0 * (DG_ORDER + 1) * (DG_ORDER + 2) * fscale[0];
  tau[_faceNum] = gtau;
}