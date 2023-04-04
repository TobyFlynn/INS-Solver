inline void fpmf_3d_calc_tau_faces(const int **order, const int *faceNum,
                              const int *fmaskR_corrected, const DG_FP *fscale,
                              const DG_FP **factor, DG_FP *tau) {
  const int p = order[0][0];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

  const int *fmask  = &FMASK[(p - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * dg_npf];

  DG_FP gtau = 0.0;
  for(int j = 0; j < dg_npf; j++) {
    DG_FP tmp = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(fscale[0] * factor[0][fmaskL[j]], fscale[1] * factor[1][fmaskR_corrected[j]]);
    gtau = fmax(gtau, tmp);
  }
  *tau = gtau;
}
