inline void ins_2d_advec_sc_rhs_oi_1(const int *faceNum, const bool *reverse,
                    const DG_FP *nx, const DG_FP *ny, const DG_FP *fscale,
                    const DG_FP **us, const DG_FP **vs, const DG_FP **ub,
                    const DG_FP **vb, DG_FP **fU, DG_FP **fV) {
  const bool rev = *reverse;
  const int edgeL = faceNum[0];
  const int edgeR = faceNum[1];
  const int *fmaskL = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeL * DG_NPF];
  const int *fmaskR = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeR * DG_NPF];
  const int fIndL = edgeL * DG_NPF;
  const int fIndR = edgeR * DG_NPF;
  const int fIndCubL = edgeL * DG_CUB_SURF_2D_NP;
  const int fIndCubR = edgeR * DG_CUB_SURF_2D_NP;

  DG_FP mUs[DG_CUB_SURF_2D_NP], mVs[DG_CUB_SURF_2D_NP];
  DG_FP mUb[DG_CUB_SURF_2D_NP], mVb[DG_CUB_SURF_2D_NP];
  DG_FP pUs[DG_CUB_SURF_2D_NP], pVs[DG_CUB_SURF_2D_NP];
  DG_FP pUb[DG_CUB_SURF_2D_NP], pVb[DG_CUB_SURF_2D_NP];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    mUs[i] = 0.0; mVs[i] = 0.0;
    mUb[i] = 0.0; mVb[i] = 0.0;
    pUs[i] = 0.0; pVs[i] = 0.0;
    pUb[i] = 0.0; pVb[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];

    DG_FP _usL = us[0][fmaskL_ind];
    DG_FP _vsL = vs[0][fmaskL_ind];
    DG_FP _ubL = ub[0][fmaskL_ind];
    DG_FP _vbL = vb[0][fmaskL_ind];
    DG_FP _usR = us[1][fmaskR_ind];
    DG_FP _vsR = vs[1][fmaskR_ind];
    DG_FP _ubR = ub[1][fmaskR_ind];
    DG_FP _vbR = vb[1][fmaskR_ind];

    for(int j = 0; j < DG_CUB_SURF_2D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCubL + j, fIndL + i, DG_CUB_SURF_2D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf2d_Interp_kernel[ind];
      mUs[j] += mat_val * _usL;
      mVs[j] += mat_val * _vsL;
      mUb[j] += mat_val * _ubL;
      mVb[j] += mat_val * _vbL;
      pUs[j] += mat_val * _usR;
      pVs[j] += mat_val * _vsR;
      pUb[j] += mat_val * _ubR;
      pVb[j] += mat_val * _vbR;
    }
  }

  const DG_FP _nxL = nx[0];
  const DG_FP _nyL = ny[0];
  const DG_FP _fscaleL = fscale[0];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    const DG_FP _mUs = mUs[i];
    const DG_FP _mVs = mVs[i];
    const DG_FP _mUb = mUb[i];
    const DG_FP _mVb = mVb[i];
    const DG_FP _pUs = pUs[i];
    const DG_FP _pVs = pVs[i];
    const DG_FP _pUb = pUb[i];
    const DG_FP _pVb = pVb[i];

    const DG_FP velL = _nxL * _mUb + _nyL * _mVb;
    const DG_FP velR = _nxL * _pUb + _nyL * _pVb;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[0][fIndCubL + i] = 0.5 * _fscaleL * (_nxL * (_mUb * _mUs + _pUb * _pUs) + _nyL * (_mUb * _mVs + _pUb * _pVs)
          + maxvel * (_mUs - _pUs));
    fV[0][fIndCubL + i] = 0.5 * _fscaleL * (_nxL * (_mVb * _mUs + _pVb * _pUs) + _nyL * (_mVb * _mVs + _pVb * _pVs)
          + maxvel * (_mVs - _pVs));
  }

  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    mUs[i] = 0.0; mVs[i] = 0.0;
    mUb[i] = 0.0; mVb[i] = 0.0;
    pUs[i] = 0.0; pVs[i] = 0.0;
    pUb[i] = 0.0; pVb[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskR_ind = fmaskR[i];
    const int fmaskL_ind = rev ? fmaskL[DG_NPF - i - 1] : fmaskL[i];

    DG_FP _usL = us[0][fmaskL_ind];
    DG_FP _vsL = vs[0][fmaskL_ind];
    DG_FP _ubL = ub[0][fmaskL_ind];
    DG_FP _vbL = vb[0][fmaskL_ind];
    DG_FP _usR = us[1][fmaskR_ind];
    DG_FP _vsR = vs[1][fmaskR_ind];
    DG_FP _ubR = ub[1][fmaskR_ind];
    DG_FP _vbR = vb[1][fmaskR_ind];

    for(int j = 0; j < DG_CUB_SURF_2D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCubR + j, fIndR + i, DG_CUB_SURF_2D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf2d_Interp_kernel[ind];
      pUs[j] += mat_val * _usL;
      pVs[j] += mat_val * _vsL;
      pUb[j] += mat_val * _ubL;
      pVb[j] += mat_val * _vbL;
      mUs[j] += mat_val * _usR;
      mVs[j] += mat_val * _vsR;
      mUb[j] += mat_val * _ubR;
      mVb[j] += mat_val * _vbR;
    }
  }

  const DG_FP _nxR = nx[1];
  const DG_FP _nyR = ny[1];
  const DG_FP _fscaleR = fscale[1];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    const DG_FP _mUs = mUs[i];
    const DG_FP _mVs = mVs[i];
    const DG_FP _mUb = mUb[i];
    const DG_FP _mVb = mVb[i];
    const DG_FP _pUs = pUs[i];
    const DG_FP _pVs = pVs[i];
    const DG_FP _pUb = pUb[i];
    const DG_FP _pVb = pVb[i];

    const DG_FP velL = _nxR * _mUb + _nyR * _mVb;
    const DG_FP velR = _nxR * _pUb + _nyR * _pVb;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[1][fIndCubR + i] = 0.5 * _fscaleR * (_nxR * (_mUb * _mUs + _pUb * _pUs) + _nyR * (_mUb * _mVs + _pUb * _pVs)
          + maxvel * (_mUs - _pUs));
    fV[1][fIndCubR + i] = 0.5 * _fscaleR * (_nxR * (_mVb * _mUs + _pVb * _pUs) + _nyR * (_mVb * _mVs + _pVb * _pVs)
          + maxvel * (_mVs - _pVs));
  }
}
