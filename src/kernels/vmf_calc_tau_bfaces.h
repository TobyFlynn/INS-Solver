inline void vmf_calc_tau_bfaces(const int *faceNum, const DG_FP *fscale, DG_FP *tau) {
  const DG_FP gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 2) * fscale[0];
  const int _faceNum = faceNum[0];
  tau[_faceNum] = gtau;
}