inline void fpmf_3d_calc_tau_faces(const int **p, const int *faceNum,
                             const int *fmaskR_corrected, const DG_FP *fscale,
                             const DG_FP **factor, DG_FP *tau) {
  const int dg_np  = DG_CONSTANTS[(p[0][0] - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(p[0][0] - 1) * DG_NUM_CONSTANTS + 1];
  const int *fmask  = &FMASK[(p[0][0] - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * dg_npf];

  const DG_FP tau_order = (DG_FP) p[0][0]; // (DG_FP) DG_ORDER;
  DG_FP tmp_tau = 0.0;
  for(int j = 0; j < dg_npf; j++) {
    const int fmaskL_ind = fmaskL[j];
    const int fmaskR_ind = fmaskR_corrected[j];
    DG_FP tmp = 2.0 * (tau_order + 1) * (tau_order + 2) * fmax(fscale[0] * factor[0][fmaskL_ind], fscale[1] * factor[1][fmaskR_ind]);
    tmp_tau = fmax(tmp_tau, tmp);
  }

  *tau = tmp_tau;
}
