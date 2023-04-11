inline void fpmf_3d_calc_tau_bflux(const int *p, const int *faceL,
                          const int *faceNums, const int *fmaskL_corrected,
                          const int *fmaskR_corrected, const DG_FP *fscale,
                          const DG_FP **factor, DG_FP *tau) {
  const int dg_np  = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 1];
  const int *fmask  = &FMASK[(*p - 1) * 4 * DG_NPF];
  const int faceInd = *faceL ? 0 : 1;

  const DG_FP tau_order = (DG_FP) *p; // (DG_FP) DG_ORDER;

  const int findL = faceNums[faceInd] * dg_npf;
  const int *fmaskL = &fmask[faceNums[faceInd] * dg_npf];
  const int *fmaskR = *faceL ? fmaskR_corrected : fmaskL_corrected;

  DG_FP gtau = 0.0;
  for(int j = 0; j < dg_npf; j++) {
    DG_FP tmp = 2.0 * (tau_order + 1) * (tau_order + 2) * fmax(fscale[0] * factor[0][fmaskL[j]], fscale[1] * factor[1][fmaskR[j]]);
    gtau = fmax(gtau, tmp);
  }

  *tau = gtau;
}
