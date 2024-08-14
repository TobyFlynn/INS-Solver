inline void euler_2d_surf_oi(const int *faceNum, const bool *reverse, const DG_FP *nx, 
                             const DG_FP *ny, const DG_FP *fscale, const DG_FP **q0, 
                             const DG_FP **q1, const DG_FP **q2, const DG_FP **q3, 
                             DG_FP **out0, DG_FP **out1, DG_FP **out2, DG_FP **out3) {
  const bool rev = *reverse;
  const int edgeL = faceNum[0];
  const int edgeR = faceNum[1];
  const int *fmaskL = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeL * DG_NPF];
  const int *fmaskR = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeR * DG_NPF];
  const int fIndL = edgeL * DG_NPF;
  const int fIndR = edgeR * DG_NPF;
  const int fIndCubL = edgeL * DG_CUB_SURF_2D_NP;
  const int fIndCubR = edgeR * DG_CUB_SURF_2D_NP;

  DG_FP mQ0[DG_CUB_SURF_2D_NP], mQ1[DG_CUB_SURF_2D_NP];
  DG_FP mQ2[DG_CUB_SURF_2D_NP], mQ3[DG_CUB_SURF_2D_NP];
  DG_FP pQ0[DG_CUB_SURF_2D_NP], pQ1[DG_CUB_SURF_2D_NP];
  DG_FP pQ2[DG_CUB_SURF_2D_NP], pQ3[DG_CUB_SURF_2D_NP];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    mQ0[i] = 0.0; mQ1[i] = 0.0;
    mQ2[i] = 0.0; mQ3[i] = 0.0;
    pQ0[i] = 0.0; pQ1[i] = 0.0;
    pQ2[i] = 0.0; pQ3[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];

    DG_FP _q0L = q0[0][fmaskL_ind];
    DG_FP _q1L = q1[0][fmaskL_ind];
    DG_FP _q2L = q2[0][fmaskL_ind];
    DG_FP _q3L = q3[0][fmaskL_ind];
    DG_FP _q0R = q0[1][fmaskR_ind];
    DG_FP _q1R = q1[1][fmaskR_ind];
    DG_FP _q2R = q2[1][fmaskR_ind];
    DG_FP _q3R = q3[1][fmaskR_ind];

    for(int j = 0; j < DG_CUB_SURF_2D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCubL + j, fIndL + i, DG_CUB_SURF_2D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf2d_Interp_kernel[ind];
      mQ0[j] += mat_val * _q0L;
      mQ1[j] += mat_val * _q1L;
      mQ2[j] += mat_val * _q2L;
      mQ3[j] += mat_val * _q3L;
      pQ0[j] += mat_val * _q0R;
      pQ1[j] += mat_val * _q1R;
      pQ2[j] += mat_val * _q2R;
      pQ3[j] += mat_val * _q3R;
    }
  }

  const DG_FP _nxL = nx[0];
  const DG_FP _nyL = ny[0];
  const DG_FP _fscaleL = fscale[0];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    const DG_FP mRho = mQ0[i];
    const DG_FP mRhou = mQ1[i];
    const DG_FP mRhov = mQ2[i];
    const DG_FP mEnergy = mQ3[i];
    const DG_FP mU = mRhou / mRho;
    const DG_FP mV = mRhov / mRho;
    const DG_FP mP = (1.4 - 1.0) * (mEnergy - 0.5 * (mRhou * mU + mRhov * mV));

    const DG_FP pRho = pQ0[i];
    const DG_FP pRhou = pQ1[i];
    const DG_FP pRhov = pQ2[i];
    const DG_FP pEnergy = pQ3[i];
    const DG_FP pU = pRhou / pRho;
    const DG_FP pV = pRhov / pRho;
    const DG_FP pP = (1.4 - 1.0) * (pEnergy - 0.5 * (pRhou * pU + pRhov * pV));

    const DG_FP lambda = fmax(sqrt(mU * mU + mV * mV) + sqrt(fabs(1.4 * mP / mRho)), sqrt(pU * pU + pV * pV) + sqrt(fabs(1.4 * pP / pRho)));

    out0[0][fIndCubL + i] = 0.5 * _fscaleL * (_nxL * (mRhou + pRhou) + _nyL * (mRhov + pRhov) + lambda * (mRho - pRho));
    out1[0][fIndCubL + i] = 0.5 * _fscaleL * (_nxL * (mRhou * mU + mP + pRhou * pU + pP) + _nyL * (mRhou * mV + pRhou * pV) + lambda * (mRhou - pRhou));
    out2[0][fIndCubL + i] = 0.5 * _fscaleL * (_nxL * (mRhou * mV + pRhou * pV) + _nyL * (mRhov * mV + mP + pRhov * pV + pP) + lambda * (mRhov - pRhov));
    out3[0][fIndCubL + i] = 0.5 * _fscaleL * (_nxL * (mU * (mEnergy + mP) + pU * (pEnergy + pP)) + _nyL * (mV * (mEnergy + mP) + pV * (pEnergy + pP)) + lambda * (mEnergy - pEnergy));
  }

  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    mQ0[i] = 0.0; mQ1[i] = 0.0;
    mQ2[i] = 0.0; mQ3[i] = 0.0;
    pQ0[i] = 0.0; pQ1[i] = 0.0;
    pQ2[i] = 0.0; pQ3[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskR_ind = fmaskR[i];
    const int fmaskL_ind = rev ? fmaskL[DG_NPF - i - 1] : fmaskL[i];

    DG_FP _q0L = q0[0][fmaskL_ind];
    DG_FP _q1L = q1[0][fmaskL_ind];
    DG_FP _q2L = q2[0][fmaskL_ind];
    DG_FP _q3L = q3[0][fmaskL_ind];
    DG_FP _q0R = q0[1][fmaskR_ind];
    DG_FP _q1R = q1[1][fmaskR_ind];
    DG_FP _q2R = q2[1][fmaskR_ind];
    DG_FP _q3R = q3[1][fmaskR_ind];

    for(int j = 0; j < DG_CUB_SURF_2D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCubR + j, fIndR + i, DG_CUB_SURF_2D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf2d_Interp_kernel[ind];
      mQ0[j] += mat_val * _q0R;
      mQ1[j] += mat_val * _q1R;
      mQ2[j] += mat_val * _q2R;
      mQ3[j] += mat_val * _q3R;
      pQ0[j] += mat_val * _q0L;
      pQ1[j] += mat_val * _q1L;
      pQ2[j] += mat_val * _q2L;
      pQ3[j] += mat_val * _q3L;
    }
  }

  const DG_FP _nxR = nx[1];
  const DG_FP _nyR = ny[1];
  const DG_FP _fscaleR = fscale[1];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    const DG_FP mRho = mQ0[i];
    const DG_FP mRhou = mQ1[i];
    const DG_FP mRhov = mQ2[i];
    const DG_FP mEnergy = mQ3[i];
    const DG_FP mU = mRhou / mRho;
    const DG_FP mV = mRhov / mRho;
    const DG_FP mP = (1.4 - 1.0) * (mEnergy - 0.5 * (mRhou * mU + mRhov * mV));

    const DG_FP pRho = pQ0[i];
    const DG_FP pRhou = pQ1[i];
    const DG_FP pRhov = pQ2[i];
    const DG_FP pEnergy = pQ3[i];
    const DG_FP pU = pRhou / pRho;
    const DG_FP pV = pRhov / pRho;
    const DG_FP pP = (1.4 - 1.0) * (pEnergy - 0.5 * (pRhou * pU + pRhov * pV));

    const DG_FP lambda = fmax(sqrt(mU * mU + mV * mV) + sqrt(fabs(1.4 * mP / mRho)), sqrt(pU * pU + pV * pV) + sqrt(fabs(1.4 * pP / pRho)));

    out0[1][fIndCubR + i] = 0.5 * _fscaleR * (_nxR * (mRhou + pRhou) + _nyR * (mRhov + pRhov) + lambda * (mRho - pRho));
    out1[1][fIndCubR + i] = 0.5 * _fscaleR * (_nxR * (mRhou * mU + mP + pRhou * pU + pP) + _nyR * (mRhou * mV + pRhou * pV) + lambda * (mRhou - pRhou));
    out2[1][fIndCubR + i] = 0.5 * _fscaleR * (_nxR * (mRhou * mV + pRhou * pV) + _nyR * (mRhov * mV + mP + pRhov * pV + pP) + lambda * (mRhov - pRhov));
    out3[1][fIndCubR + i] = 0.5 * _fscaleR * (_nxR * (mU * (mEnergy + mP) + pU * (pEnergy + pP)) + _nyR * (mV * (mEnergy + mP) + pV * (pEnergy + pP)) + lambda * (mEnergy - pEnergy));
  }
}
