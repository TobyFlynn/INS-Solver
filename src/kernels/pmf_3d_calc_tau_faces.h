inline void pmf_3d_calc_tau_faces(const int **order, const int *faceNum,
                                  const DG_FP *fscale, DG_FP **tau) {
  const int p = order[0][0];
  const DG_FP gtau = 2.0 * (p + 1) * (p + 2) * fmax(fscale[0], fscale[1]);
  const int faceNumL = faceNum[0];
  const int faceNumR = faceNum[1];
  tau[0][faceNumL] = gtau;
  tau[1][faceNumR] = gtau;
}
