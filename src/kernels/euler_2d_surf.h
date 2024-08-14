inline void euler_2d_surf(const int *faceNum, const bool *reverse, const DG_FP *nx, 
                          const DG_FP *ny, const DG_FP *fscale, const DG_FP **q0, 
                          const DG_FP **q1, const DG_FP **q2, const DG_FP **q3, 
                          DG_FP **out0, DG_FP **out1, DG_FP **out2, DG_FP **out3) {
  const bool rev = *reverse;
  const int edgeL = faceNum[0];
  const int edgeR = faceNum[1];
  const int *fmaskL = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeL * DG_NPF];
  const int *fmaskR = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeR * DG_NPF];

  const DG_FP _nxL = nx[0];
  const DG_FP _nyL = ny[0];
  const DG_FP _fscaleL = fscale[0];
  const DG_FP _nxR = nx[1];
  const DG_FP _nyR = ny[1];
  const DG_FP _fscaleR = fscale[1];
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];

    const int fIndL = edgeL * DG_NPF + i;
    const int fIndR = rev ? edgeR * DG_NPF + DG_NPF - i - 1 : edgeR * DG_NPF + i;

    const DG_FP mRho = q0[0][fmaskL_ind];
    const DG_FP mRhou = q1[0][fmaskL_ind];
    const DG_FP mRhov = q2[0][fmaskL_ind];
    const DG_FP mEnergy = q3[0][fmaskL_ind];
    const DG_FP mU = mRhou / mRho;
    const DG_FP mV = mRhov / mRho;
    const DG_FP mP = (1.4 - 1.0) * (mEnergy - 0.5 * (mRhou * mU + mRhov * mV));

    const DG_FP pRho = q0[1][fmaskR_ind];
    const DG_FP pRhou = q1[1][fmaskR_ind];
    const DG_FP pRhov = q2[1][fmaskR_ind];
    const DG_FP pEnergy = q3[1][fmaskR_ind];
    const DG_FP pU = pRhou / pRho;
    const DG_FP pV = pRhov / pRho;
    const DG_FP pP = (1.4 - 1.0) * (pEnergy - 0.5 * (pRhou * pU + pRhov * pV));

    const DG_FP lambda = fmax(sqrt(mU * mU + mV * mV) + sqrt(fabs(1.4 * mP / mRho)), sqrt(pU * pU + pV * pV) + sqrt(fabs(1.4 * pP / pRho)));

    out0[0][fIndL] = 0.5 * _fscaleL * (_nxL * (mRhou + pRhou) + _nyL * (mRhov + pRhov) + lambda * (mRho - pRho));
    out1[0][fIndL] = 0.5 * _fscaleL * (_nxL * (mRhou * mU + mP + pRhou * pU + pP) + _nyL * (mRhou * mV + pRhou * pV) + lambda * (mRhou - pRhou));
    out2[0][fIndL] = 0.5 * _fscaleL * (_nxL * (mRhou * mV + pRhou * pV) + _nyL * (mRhov * mV + mP + pRhov * pV + pP) + lambda * (mRhov - pRhov));
    out3[0][fIndL] = 0.5 * _fscaleL * (_nxL * (mU * (mEnergy + mP) + pU * (pEnergy + pP)) + _nyL * (mV * (mEnergy + mP) + pV * (pEnergy + pP)) + lambda * (mEnergy - pEnergy));

    out0[1][fIndR] = 0.5 * _fscaleR * (_nxR * (mRhou + pRhou) + _nyR * (mRhov + pRhov) + lambda * (pRho - mRho));
    out1[1][fIndR] = 0.5 * _fscaleR * (_nxR * (mRhou * mU + mP + pRhou * pU + pP) + _nyR * (mRhou * mV + pRhou * pV) + lambda * (pRhou - mRhou));
    out2[1][fIndR] = 0.5 * _fscaleR * (_nxR * (mRhou * mV + pRhou * pV) + _nyR * (mRhov * mV + mP + pRhov * pV + pP) + lambda * (pRhov - mRhov));
    out3[1][fIndR] = 0.5 * _fscaleR * (_nxR * (mU * (mEnergy + mP) + pU * (pEnergy + pP)) + _nyR * (mV * (mEnergy + mP) + pV * (pEnergy + pP)) + lambda * (pEnergy - mEnergy));
  }
}
