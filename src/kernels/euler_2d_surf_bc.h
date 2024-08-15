inline void euler_2d_surf_bc(const int *edgeNum, const DG_FP *nx,
                const DG_FP *ny, const DG_FP *fscale, const DG_FP *q0,
                const DG_FP *q1, const DG_FP *q2, const DG_FP *q3,
                DG_FP *out0, DG_FP *out1, DG_FP *out2, DG_FP *out3) {
  // Work out which edge for each element
  const int edge = edgeNum[0];
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];
  const DG_FP _nx = nx[0];
  const DG_FP _ny = ny[0];
  const DG_FP _fscale = fscale[0];
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];

    const DG_FP mRho = q0[fmask_ind];
    const DG_FP mRhou = q1[fmask_ind];
    const DG_FP mRhov = q2[fmask_ind];
    const DG_FP mEnergy = q3[fmask_ind];
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

    out0[edge * DG_NPF + i] += 0.5 * _fscale * (_nx * (mRhou + pRhou) + _ny * (mRhov + pRhov) + lambda * (mRho - pRho));
    out1[edge * DG_NPF + i] += 0.5 * _fscale * (_nx * (mRhou * mU + mP + pRhou * pU + pP) + _ny * (mRhou * mV + pRhou * pV) + lambda * (mRhou - pRhou));
    out2[edge * DG_NPF + i] += 0.5 * _fscale * (_nx * (mRhou * mV + pRhou * pV) + _ny * (mRhov * mV + mP + pRhov * pV + pP) + lambda * (mRhov - pRhov));
    out3[edge * DG_NPF + i] += 0.5 * _fscale * (_nx * (mU * (mEnergy + mP) + pU * (pEnergy + pP)) + _ny * (mV * (mEnergy + mP) + pV * (pEnergy + pP)) + lambda * (mEnergy - pEnergy));
  }
}