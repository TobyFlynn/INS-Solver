inline void pmf_3d_calc_tau_bfaces_sp(const int *order, const int *faceNum,
                                   const DG_FP *fscale, float *tau) {
  const int p = order[0];
  const DG_FP gtau = 2.0 * (p + 1) * (p + 2) * fscale[0];
  const int _faceNum = faceNum[0];
  tau[_faceNum] = (float)gtau;
}
