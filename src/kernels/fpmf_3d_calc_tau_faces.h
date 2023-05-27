inline void fpmf_3d_calc_tau_faces(const int **order, const int *faceNum,
                                   const int *fmaskR_corrected,
                                   const DG_FP *fscale, const DG_FP **factor,
                                   DG_FP **tau) {
  const int p = order[0][0];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];
  DG_FP gtau = 0.0;
  const int faceNumL = faceNum[0];
  for(int j = 0; j < dg_npf; j++) {
    const int fmaskL_ind = FMASK[(p - 1) * 4 * DG_NPF + faceNumL * dg_npf + j];
    const int fmaskR_ind = fmaskR_corrected[j];
    DG_FP tmp = 2.0 * (p + 1) * (p + 2) * fmax(fscale[0] * factor[0][fmaskL_ind], fscale[1] * factor[1][fmaskR_ind]);
    gtau = fmax(gtau, tmp);
  }
  const int faceNumR = faceNum[1];
  tau[0][faceNumL] = gtau;
  tau[1][faceNumR] = gtau;
}
