inline void fvmf_3d_calc_tau_faces(const int *faceNum, const int *fmaskR_corrected,
                                   const DG_FP *fscale, const DG_FP **factor, DG_FP **tau) {
  DG_FP gtau = 0.0;
  const int faceNumL = faceNum[0];
  for(int j = 0; j < DG_NPF; j++) {
    const int fmaskL_ind = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + faceNumL * DG_NPF + j];
    const int fmaskR_ind = fmaskR_corrected[j];
    DG_FP tmp = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(fscale[0] * factor[0][fmaskL_ind], fscale[1] * factor[1][fmaskR_ind]);
    gtau = fmax(gtau, tmp);
  }
  const int faceNumR = faceNum[1];
  tau[0][faceNumL] = gtau;
  tau[1][faceNumR] = gtau;
}