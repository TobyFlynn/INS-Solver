inline void euler_2d_surf_bc_oi(const int *edgeNum, const DG_FP *nx,
                const DG_FP *ny, const DG_FP *fscale, const DG_FP *q0,
                const DG_FP *q1, const DG_FP *q2, const DG_FP *q3,
                DG_FP *out0, DG_FP *out1, DG_FP *out2, DG_FP *out3) {
  const int edge = edgeNum[0];
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];
  const DG_FP _nx = nx[0];
  const DG_FP _ny = ny[0];
  const DG_FP _fscale = fscale[0];
  const int fIndCub = edge * DG_CUB_SURF_2D_NP;
  const int fInd = edge * DG_NPF;

  DG_FP mQ0[DG_CUB_SURF_2D_NP], mQ1[DG_CUB_SURF_2D_NP];
  DG_FP mQ2[DG_CUB_SURF_2D_NP], mQ3[DG_CUB_SURF_2D_NP];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    mQ0[i] = 0.0; mQ1[i] = 0.0;
    mQ2[i] = 0.0; mQ3[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];

    DG_FP _q0L = q0[fmask_ind];
    DG_FP _q1L = q1[fmask_ind];
    DG_FP _q2L = q2[fmask_ind];
    DG_FP _q3L = q3[fmask_ind];

    for(int j = 0; j < DG_CUB_SURF_2D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCub + j, fInd + i, DG_CUB_SURF_2D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf2d_Interp_kernel[ind];
      mQ0[j] += mat_val * _q0L;
      mQ1[j] += mat_val * _q1L;
      mQ2[j] += mat_val * _q2L;
      mQ3[j] += mat_val * _q3L;
    }
  }

  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    const DG_FP mRho = mQ0[i];
    const DG_FP mRhou = mQ1[i];
    const DG_FP mRhov = mQ2[i];
    const DG_FP mEnergy = mQ3[i];
    const DG_FP mU = mRhou / mRho;
    const DG_FP mV = mRhov / mRho;
    const DG_FP mP = (1.4 - 1.0) * (mEnergy - 0.5 * (mRhou * mU + mRhov * mV));

    const DG_FP pRho = 1.0;
    const DG_FP pRhou = 0.0;
    const DG_FP pRhov = 1.0;
    const DG_FP pP = 1.0 / (1.4 * 0.4 * 0.4);
    const DG_FP pU = pRhou / pRho;
    const DG_FP pV = pRhov / pRho;
    const DG_FP pEnergy = (pP / 0.4) + 0.5 * pRho * (pU * pU + pV * pV);

    const DG_FP lambda = fmax(sqrt(mU * mU + mV * mV) + sqrt(fabs(1.4 * mP / mRho)), sqrt(pU * pU + pV * pV) + sqrt(fabs(1.4 * pP / pRho)));

    out0[fIndCub + i] = 0.5 * _fscale * (_nx * (mRhou + pRhou) + _ny * (mRhov + pRhov) + lambda * (mRho - pRho));
    out1[fIndCub + i] = 0.5 * _fscale * (_nx * (mRhou * mU + mP + pRhou * pU + pP) + _ny * (mRhou * mV + pRhou * pV) + lambda * (mRhou - pRhou));
    out2[fIndCub + i] = 0.5 * _fscale * (_nx * (mRhou * mV + pRhou * pV) + _ny * (mRhov * mV + mP + pRhov * pV + pP) + lambda * (mRhov - pRhov));
    out3[fIndCub + i] = 0.5 * _fscale * (_nx * (mU * (mEnergy + mP) + pU * (pEnergy + pP)) + _ny * (mV * (mEnergy + mP) + pV * (pEnergy + pP)) + lambda * (mEnergy - pEnergy));
  }
}