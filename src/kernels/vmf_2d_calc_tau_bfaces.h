inline void vmf_2d_calc_tau_bfaces(const int *order, const int *faceNum,
                                   const DG_FP *fscale, DG_FP *tau) {
  const int p = order[0];
  const DG_FP gtau = 2.0 * (p + 1) * (p + 2) * fscale[0];
  const int _faceNum = faceNum[0];
  tau[_faceNum] = gtau;
}