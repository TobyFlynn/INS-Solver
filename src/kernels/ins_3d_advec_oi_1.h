inline void ins_3d_advec_oi_1(const int *faceNum, const int *fmaskL_corrected,
                        const int *fmaskR_corrected, const DG_FP *nx,
                        const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                        const DG_FP **u, const DG_FP **v, const DG_FP **w,
                        DG_FP **fU, DG_FP **fV, DG_FP **fW) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];
  const int fIndL = faceNum[0] * DG_NPF;
  const int fIndR = faceNum[1] * DG_NPF;
  const int fIndCubL = faceNum[0] * DG_CUB_SURF_3D_NP;
  const int fIndCubR = faceNum[1] * DG_CUB_SURF_3D_NP;

  DG_FP mU[DG_CUB_SURF_3D_NP], mV[DG_CUB_SURF_3D_NP], mW[DG_CUB_SURF_3D_NP];
  DG_FP pU[DG_CUB_SURF_3D_NP], pV[DG_CUB_SURF_3D_NP], pW[DG_CUB_SURF_3D_NP];
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    mU[i] = 0.0; mV[i] = 0.0; mW[i] = 0.0;
    pU[i] = 0.0; pV[i] = 0.0; pW[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];

    DG_FP _uL = u[0][fmaskL_ind];
    DG_FP _vL = v[0][fmaskL_ind];
    DG_FP _wL = w[0][fmaskL_ind];
    DG_FP _uR = u[1][fmaskR_ind];
    DG_FP _vR = v[1][fmaskR_ind];
    DG_FP _wR = w[1][fmaskR_ind];

    for(int j = 0; j < DG_CUB_SURF_3D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCubL + j, fIndL + i, DG_CUB_SURF_3D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf3d_Interp_kernel[ind];
      mU[j] += mat_val * _uL;
      mV[j] += mat_val * _vL;
      mW[j] += mat_val * _wL;
      pU[j] += mat_val * _uR;
      pV[j] += mat_val * _vR;
      pW[j] += mat_val * _wR;
    }
  }

  const DG_FP _nxL = nx[0];
  const DG_FP _nyL = ny[0];
  const DG_FP _nzL = nz[0];
  const DG_FP _fscaleL = fscale[0];
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    const DG_FP _mU = mU[i];
    const DG_FP _mV = mV[i];
    const DG_FP _mW = mW[i];
    const DG_FP _pU = pU[i];
    const DG_FP _pV = pV[i];
    const DG_FP _pW = pW[i];

    const DG_FP velL = _nxL * _mU + _nyL * _mV + _nzL * _mW;
    const DG_FP velR = _nxL * _pU + _nyL * _pV + _nzL * _pW;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[0][fIndCubL + i] = 0.5 * _fscaleL * (_nxL * (_mU * _mU + _pU * _pU) + _nyL * (_mU * _mV + _pU * _pV)
          + _nzL * (_mU * _mW + _pU * _pW) + maxvel * (_mU - _pU));
    fV[0][fIndCubL + i] = 0.5 * _fscaleL * (_nxL * (_mV * _mU + _pV * _pU) + _nyL * (_mV * _mV + _pV * _pV)
          + _nzL * (_mV * _mW + _pV * _pW) + maxvel * (_mV - _pV));
    fW[0][fIndCubL + i] = 0.5 * _fscaleL * (_nxL * (_mW * _mU + _pW * _pU) + _nyL * (_mW * _mV + _pW * _pV)
          + _nzL * (_mW * _mW + _pW * _pW) + maxvel * (_mW - _pW));
  }

  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    mU[i] = 0.0; mV[i] = 0.0; mW[i] = 0.0;
    pU[i] = 0.0; pV[i] = 0.0; pW[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL_corrected[i];
    const int fmaskR_ind = fmaskR[i];

    DG_FP _uL = u[0][fmaskL_ind];
    DG_FP _vL = v[0][fmaskL_ind];
    DG_FP _wL = w[0][fmaskL_ind];
    DG_FP _uR = u[1][fmaskR_ind];
    DG_FP _vR = v[1][fmaskR_ind];
    DG_FP _wR = w[1][fmaskR_ind];

    for(int j = 0; j < DG_CUB_SURF_3D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCubR + j, fIndR + i, DG_CUB_SURF_3D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf3d_Interp_kernel[ind];
      mU[j] += mat_val * _uR;
      mV[j] += mat_val * _vR;
      mW[j] += mat_val * _wR;
      pU[j] += mat_val * _uL;
      pV[j] += mat_val * _vL;
      pW[j] += mat_val * _wL;
    }
  }

  const DG_FP _nxR = nx[1];
  const DG_FP _nyR = ny[1];
  const DG_FP _nzR = nz[1];
  const DG_FP _fscaleR = fscale[1];
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    const DG_FP _mU = mU[i];
    const DG_FP _mV = mV[i];
    const DG_FP _mW = mW[i];
    const DG_FP _pU = pU[i];
    const DG_FP _pV = pV[i];
    const DG_FP _pW = pW[i];

    const DG_FP velL = _nxR * _mU + _nyR * _mV + _nzR * _mW;
    const DG_FP velR = _nxR * _pU + _nyR * _pV + _nzR * _pW;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[1][fIndCubR + i] = 0.5 * _fscaleR * (_nxR * (_mU * _mU + _pU * _pU) + _nyR * (_mU * _mV + _pU * _pV)
          + _nzR * (_mU * _mW + _pU * _pW) + maxvel * (_mU - _pU));
    fV[1][fIndCubR + i] = 0.5 * _fscaleR * (_nxR * (_mV * _mU + _pV * _pU) + _nyR * (_mV * _mV + _pV * _pV)
          + _nzR * (_mV * _mW + _pV * _pW) + maxvel * (_mV - _pV));
    fW[1][fIndCubR + i] = 0.5 * _fscaleR * (_nxR * (_mW * _mU + _pW * _pU) + _nyR * (_mW * _mV + _pW * _pV)
          + _nzR * (_mW * _mW + _pW * _pW) + maxvel * (_mW - _pW));
  }
}
