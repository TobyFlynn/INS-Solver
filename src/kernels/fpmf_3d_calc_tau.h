inline void fpmf_3d_calc_tau(const int *p, const int *faceNums,
                             const int *fmaskF, const DG_FP *fscale,
                             const DG_FP *factor, const DG_FP **factor_p,
                             DG_FP *tau) {
  const int dg_np  = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 1];
  const int *fmask  = &FMASK[(*p - 1) * 4 * DG_NPF];

  for(int i = 0; i < DG_NUM_FACES; i++) {
    const int *fmaskL = &fmask[faceNums[2 * i] * dg_npf];
    const int *fmaskR = &fmaskF[i * dg_npf];

    tau[i] = 0.0;
    for(int j = 0; j < dg_npf; j++) {
      DG_FP tmp = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(fscale[i * 2] * factor[fmaskL[j]], fscale[i * 2 + 1] * factor_p[i][fmaskR[j]]);
      tau[i] = fmax(tau[i], tmp);
    }
  }
}
