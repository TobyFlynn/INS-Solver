inline void ins_3d_advec_sc_rhs_oi_1(const int *faceNum, const int *fmaskL_corrected,
                        const int *fmaskR_corrected, const DG_FP *nx,
                        const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                        const DG_FP **us, const DG_FP **vs, const DG_FP **ws,
                        const DG_FP **ub, const DG_FP **vb, const DG_FP **wb,
                        DG_FP **fU, DG_FP **fV, DG_FP **fW) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];
  const int fIndL = faceNum[0] * DG_NPF;
  const int fIndR = faceNum[1] * DG_NPF;
  const int fIndCubL = faceNum[0] * DG_CUB_SURF_3D_NP;
  const int fIndCubR = faceNum[1] * DG_CUB_SURF_3D_NP;

  DG_FP mUs[DG_CUB_SURF_3D_NP], mVs[DG_CUB_SURF_3D_NP], mWs[DG_CUB_SURF_3D_NP];
  DG_FP mUb[DG_CUB_SURF_3D_NP], mVb[DG_CUB_SURF_3D_NP], mWb[DG_CUB_SURF_3D_NP];
  DG_FP pUs[DG_CUB_SURF_3D_NP], pVs[DG_CUB_SURF_3D_NP], pWs[DG_CUB_SURF_3D_NP];
  DG_FP pUb[DG_CUB_SURF_3D_NP], pVb[DG_CUB_SURF_3D_NP], pWb[DG_CUB_SURF_3D_NP];
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    mUs[i] = 0.0; mVs[i] = 0.0; mWs[i] = 0.0;
    mUb[i] = 0.0; mVb[i] = 0.0; mWb[i] = 0.0;
    pUs[i] = 0.0; pVs[i] = 0.0; pWs[i] = 0.0;
    pUb[i] = 0.0; pVb[i] = 0.0; pWb[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];

    DG_FP _usL = us[0][fmaskL_ind];
    DG_FP _vsL = vs[0][fmaskL_ind];
    DG_FP _wsL = ws[0][fmaskL_ind];
    DG_FP _ubL = ub[0][fmaskL_ind];
    DG_FP _vbL = vb[0][fmaskL_ind];
    DG_FP _wbL = wb[0][fmaskL_ind];
    DG_FP _usR = us[1][fmaskR_ind];
    DG_FP _vsR = vs[1][fmaskR_ind];
    DG_FP _wsR = ws[1][fmaskR_ind];
    DG_FP _ubR = ub[1][fmaskR_ind];
    DG_FP _vbR = vb[1][fmaskR_ind];
    DG_FP _wbR = wb[1][fmaskR_ind];

    for(int j = 0; j < DG_CUB_SURF_3D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCubL + j, fIndL + i, DG_CUB_SURF_3D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf3d_Interp_kernel[ind];
      mUs[j] += mat_val * _usL;
      mVs[j] += mat_val * _vsL;
      mWs[j] += mat_val * _wsL;
      mUb[j] += mat_val * _ubL;
      mVb[j] += mat_val * _vbL;
      mWb[j] += mat_val * _wbL;
      pUs[j] += mat_val * _usR;
      pVs[j] += mat_val * _vsR;
      pWs[j] += mat_val * _wsR;
      pUb[j] += mat_val * _ubR;
      pVb[j] += mat_val * _vbR;
      pWb[j] += mat_val * _wbR;
    }
  }

  const DG_FP _nxL = nx[0];
  const DG_FP _nyL = ny[0];
  const DG_FP _nzL = nz[0];
  const DG_FP _fscaleL = fscale[0];
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    const DG_FP _mUs = mUs[i];
    const DG_FP _mVs = mVs[i];
    const DG_FP _mWs = mWs[i];
    const DG_FP _pUs = pUs[i];
    const DG_FP _pVs = pVs[i];
    const DG_FP _pWs = pWs[i];
    const DG_FP _mUb = mUb[i];
    const DG_FP _mVb = mVb[i];
    const DG_FP _mWb = mWb[i];
    const DG_FP _pUb = pUb[i];
    const DG_FP _pVb = pVb[i];
    const DG_FP _pWb = pWb[i];

    // Check whether it shoudl be s instead of b here (or the max of them)
    const DG_FP velL = _nxL * _mUb + _nyL * _mVb + _nzL * _mWb;
    const DG_FP velR = _nxL * _pUb + _nyL * _pVb + _nzL * _pWb;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[0][fIndCubL + i] = 0.5 * _fscaleL * (_nxL * (_mUb * _mUs + _pUb * _pUs) + _nyL * (_mUb * _mVs + _pUb * _pVs)
          + _nzL * (_mUb * _mWs + _pUb * _pWs) + maxvel * (_mUs - _pUs));
    fV[0][fIndCubL + i] = 0.5 * _fscaleL * (_nxL * (_mVb * _mUs + _pVb * _pUs) + _nyL * (_mVb * _mVs + _pVb * _pVs)
          + _nzL * (_mVb * _mWs + _pVb * _pWs) + maxvel * (_mVs - _pVs));
    fW[0][fIndCubL + i] = 0.5 * _fscaleL * (_nxL * (_mWb * _mUs + _pWb * _pUs) + _nyL * (_mWb * _mVs + _pWb * _pVs)
          + _nzL * (_mWb * _mWs + _pWb * _pWs) + maxvel * (_mWs - _pWs));
  }

  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    mUs[i] = 0.0; mVs[i] = 0.0; mWs[i] = 0.0;
    mUb[i] = 0.0; mVb[i] = 0.0; mWb[i] = 0.0;
    pUs[i] = 0.0; pVs[i] = 0.0; pWs[i] = 0.0;
    pUb[i] = 0.0; pVb[i] = 0.0; pWb[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL_corrected[i];
    const int fmaskR_ind = fmaskR[i];

    DG_FP _usL = us[0][fmaskL_ind];
    DG_FP _vsL = vs[0][fmaskL_ind];
    DG_FP _wsL = ws[0][fmaskL_ind];
    DG_FP _ubL = ub[0][fmaskL_ind];
    DG_FP _vbL = vb[0][fmaskL_ind];
    DG_FP _wbL = wb[0][fmaskL_ind];
    DG_FP _usR = us[1][fmaskR_ind];
    DG_FP _vsR = vs[1][fmaskR_ind];
    DG_FP _wsR = ws[1][fmaskR_ind];
    DG_FP _ubR = ub[1][fmaskR_ind];
    DG_FP _vbR = vb[1][fmaskR_ind];
    DG_FP _wbR = wb[1][fmaskR_ind];

    for(int j = 0; j < DG_CUB_SURF_3D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCubR + j, fIndR + i, DG_CUB_SURF_3D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf3d_Interp_kernel[ind];
      mUs[j] += mat_val * _usR;
      mVs[j] += mat_val * _vsR;
      mWs[j] += mat_val * _wsR;
      mUb[j] += mat_val * _ubR;
      mVb[j] += mat_val * _vbR;
      mWb[j] += mat_val * _wbR;
      pUs[j] += mat_val * _usL;
      pVs[j] += mat_val * _vsL;
      pWs[j] += mat_val * _wsL;
      pUb[j] += mat_val * _ubL;
      pVb[j] += mat_val * _vbL;
      pWb[j] += mat_val * _wbL;
    }
  }

  const DG_FP _nxR = nx[1];
  const DG_FP _nyR = ny[1];
  const DG_FP _nzR = nz[1];
  const DG_FP _fscaleR = fscale[1];
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    const DG_FP _mUs = mUs[i];
    const DG_FP _mVs = mVs[i];
    const DG_FP _mWs = mWs[i];
    const DG_FP _pUs = pUs[i];
    const DG_FP _pVs = pVs[i];
    const DG_FP _pWs = pWs[i];
    const DG_FP _mUb = mUb[i];
    const DG_FP _mVb = mVb[i];
    const DG_FP _mWb = mWb[i];
    const DG_FP _pUb = pUb[i];
    const DG_FP _pVb = pVb[i];
    const DG_FP _pWb = pWb[i];

    // Check whether it shoudl be s instead of b here (or the max of them)
    const DG_FP velL = _nxR * _mUb + _nyR * _mVb + _nzR * _mWb;
    const DG_FP velR = _nxR * _pUb + _nyR * _pVb + _nzR * _pWb;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[1][fIndCubR + i] = 0.5 * _fscaleR * (_nxR * (_mUb * _mUs + _pUb * _pUs) + _nyR * (_mUb * _mVs + _pUb * _pVs)
          + _nzR * (_mUb * _mWs + _pUb * _pWs) + maxvel * (_mUs - _pUs));
    fV[1][fIndCubR + i] = 0.5 * _fscaleR * (_nxR * (_mVb * _mUs + _pVb * _pUs) + _nyR * (_mVb * _mVs + _pVb * _pVs)
          + _nzR * (_mVb * _mWs + _pVb * _pWs) + maxvel * (_mVs - _pVs));
    fW[1][fIndCubR + i] = 0.5 * _fscaleR * (_nxR * (_mWb * _mUs + _pWb * _pUs) + _nyR * (_mWb * _mVs + _pWb * _pVs)
          + _nzR * (_mWb * _mWs + _pWb * _pWs) + maxvel * (_mWs - _pWs));
  }
}
