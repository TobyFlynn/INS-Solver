inline void ins_3d_advec_sc_rhs_oi_1_new_v2(const int *faceNum, const int *fmaskL_corrected,
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

  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    DG_FP _mUs = 0.0; DG_FP _mVs = 0.0; DG_FP _mWs = 0.0;
    DG_FP _pUs = 0.0; DG_FP _pVs = 0.0; DG_FP _pWs = 0.0;
    DG_FP _mUb = 0.0; DG_FP _mVb = 0.0; DG_FP _mWb = 0.0;
    DG_FP _pUb = 0.0; DG_FP _pVb = 0.0; DG_FP _pWb = 0.0;

    for(int j = 0; j < DG_NPF; j++) {
      const int fmaskL_ind = fmaskL[i];
      const int fmaskR_ind = fmaskR_corrected[i];
      const int ind = DG_MAT_IND(fIndCubL + i, fIndL + j, DG_CUB_SURF_3D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf3d_Interp_kernel[ind];

      _mUs += mat_val * us[0][fmaskL_ind];
      _mVs += mat_val * vs[0][fmaskL_ind];
      _mWs += mat_val * ws[0][fmaskL_ind];
      _mUb += mat_val * ub[0][fmaskL_ind];
      _mVb += mat_val * vb[0][fmaskL_ind];
      _mWb += mat_val * wb[0][fmaskL_ind];
      _pUs += mat_val * us[1][fmaskR_ind];
      _pVs += mat_val * vs[1][fmaskR_ind];
      _pWs += mat_val * ws[1][fmaskR_ind];
      _pUb += mat_val * ub[1][fmaskR_ind];
      _pVb += mat_val * vb[1][fmaskR_ind];
      _pWb += mat_val * wb[1][fmaskR_ind];
    }

    const DG_FP _nx = nx[0];
    const DG_FP _ny = ny[0];
    const DG_FP _nz = nz[0];
    const DG_FP _fscale = fscale[0];

    // Check whether it shoudl be s instead of b here (or the max of them)
    const DG_FP velL = _nx * _mUb + _ny * _mVb + _nz * _mWb;
    const DG_FP velR = _nx * _pUb + _ny * _pVb + _nz * _pWb;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[0][fIndCubL + i] = 0.5 * _fscale * (_nx * (_mUb * _mUs + _pUb * _pUs) + _ny * (_mUb * _mVs + _pUb * _pVs)
          + _nz * (_mUb * _mWs + _pUb * _pWs) + maxvel * (_mUs - _pUs));
    fV[0][fIndCubL + i] = 0.5 * _fscale * (_nx * (_mVb * _mUs + _pVb * _pUs) + _ny * (_mVb * _mVs + _pVb * _pVs)
          + _nz * (_mVb * _mWs + _pVb * _pWs) + maxvel * (_mVs - _pVs));
    fW[0][fIndCubL + i] = 0.5 * _fscale * (_nx * (_mWb * _mUs + _pWb * _pUs) + _ny * (_mWb * _mVs + _pWb * _pVs)
          + _nz * (_mWb * _mWs + _pWb * _pWs) + maxvel * (_mWs - _pWs));
  }

  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    DG_FP _mUs = 0.0; DG_FP _mVs = 0.0; DG_FP _mWs = 0.0;
    DG_FP _pUs = 0.0; DG_FP _pVs = 0.0; DG_FP _pWs = 0.0;
    DG_FP _mUb = 0.0; DG_FP _mVb = 0.0; DG_FP _mWb = 0.0;
    DG_FP _pUb = 0.0; DG_FP _pVb = 0.0; DG_FP _pWb = 0.0;

    for(int j = 0; j < DG_NPF; j++) {
      const int fmaskL_ind = fmaskL_corrected[i];
      const int fmaskR_ind = fmaskR[i];
      const int ind = DG_MAT_IND(fIndCubR + i, fIndR + j, DG_CUB_SURF_3D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf3d_Interp_kernel[ind];

      _mUs += mat_val * us[1][fmaskR_ind];
      _mVs += mat_val * vs[1][fmaskR_ind];
      _mWs += mat_val * ws[1][fmaskR_ind];
      _mUb += mat_val * ub[1][fmaskR_ind];
      _mVb += mat_val * vb[1][fmaskR_ind];
      _mWb += mat_val * wb[1][fmaskR_ind];
      _pUs += mat_val * us[0][fmaskL_ind];
      _pVs += mat_val * vs[0][fmaskL_ind];
      _pWs += mat_val * ws[0][fmaskL_ind];
      _pUb += mat_val * ub[0][fmaskL_ind];
      _pVb += mat_val * vb[0][fmaskL_ind];
      _pWb += mat_val * wb[0][fmaskL_ind];
    }

    const DG_FP _nx = nx[1];
    const DG_FP _ny = ny[1];
    const DG_FP _nz = nz[1];
    const DG_FP _fscale = fscale[1];

    // Check whether it shoudl be s instead of b here (or the max of them)
    const DG_FP velL = _nx * _mUb + _ny * _mVb + _nz * _mWb;
    const DG_FP velR = _nx * _pUb + _ny * _pVb + _nz * _pWb;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[1][fIndCubR + i] = 0.5 * _fscale * (_nx * (_mUb * _mUs + _pUb * _pUs) + _ny * (_mUb * _mVs + _pUb * _pVs)
          + _nz * (_mUb * _mWs + _pUb * _pWs) + maxvel * (_mUs - _pUs));
    fV[1][fIndCubR + i] = 0.5 * _fscale * (_nx * (_mVb * _mUs + _pVb * _pUs) + _ny * (_mVb * _mVs + _pVb * _pVs)
          + _nz * (_mVb * _mWs + _pVb * _pWs) + maxvel * (_mVs - _pVs));
    fW[1][fIndCubR + i] = 0.5 * _fscale * (_nx * (_mWb * _mUs + _pWb * _pUs) + _ny * (_mWb * _mVs + _pWb * _pVs)
          + _nz * (_mWb * _mWs + _pWb * _pWs) + maxvel * (_mWs - _pWs));
  }
}
