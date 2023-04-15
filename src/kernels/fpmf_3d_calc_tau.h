inline void fpmf_3d_calc_tau(const int *p, const int *faceNums,
                             const int *fmaskF, const DG_FP *fscale,
                             const DG_FP *factor, const DG_FP **factor_p,
                             DG_FP *tau) {
  const int dg_np  = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 1];
  const int *fmask  = &FMASK[(*p - 1) * 4 * DG_NPF];

  const DG_FP tau_order = (DG_FP) *p; // (DG_FP) DG_ORDER;

  for(int i = 0; i < DG_NUM_FACES; i++) {
    const int faceNum = faceNums[2 * i];
    const int *fmaskL = &fmask[faceNum * dg_npf];

    tau[i] = 0.0;
    for(int j = 0; j < dg_npf; j++) {
      const int fmaskL_ind = fmaskL[j];
      const int fmaskR_ind = fmaskF[i * dg_npf + j];
      DG_FP tmp = 2.0 * (tau_order + 1) * (tau_order + 2) * fmax(fscale[i * 2] * factor[fmaskL_ind], fscale[i * 2 + 1] * factor_p[i][fmaskR_ind]);
      tau[i] = fmax(tau[i], tmp);
    }
  }
}
