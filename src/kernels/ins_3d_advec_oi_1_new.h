inline void ins_3d_advec_oi_1_new(const int *faceNum, const int *fmaskL_corrected,
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

  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    const DG_FP _nx = nx[0];
    const DG_FP _ny = ny[0];
    const DG_FP _nz = nz[0];
    const DG_FP _fscale = fscale[0];

    const DG_FP _mU = mU[i];
    const DG_FP _mV = mV[i];
    const DG_FP _mW = mW[i];
    const DG_FP _pU = pU[i];
    const DG_FP _pV = pV[i];
    const DG_FP _pW = pW[i];

    const DG_FP velL = _nx * _mU + _ny * _mV + _nz * _mW;
    const DG_FP velR = _nx * _pU + _ny * _pV + _nz * _pW;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[0][fIndCubL + i] = 0.5 * _fscale * (_nx * (_mU * _mU + _pU * _pU) + _ny * (_mU * _mV + _pU * _pV)
          + _nz * (_mU * _mW + _pU * _pW) + maxvel * (_mU - _pU));
    fV[0][fIndCubL + i] = 0.5 * _fscale * (_nx * (_mV * _mU + _pV * _pU) + _ny * (_mV * _mV + _pV * _pV)
          + _nz * (_mV * _mW + _pV * _pW) + maxvel * (_mV - _pV));
    fW[0][fIndCubL + i] = 0.5 * _fscale * (_nx * (_mW * _mU + _pW * _pU) + _ny * (_mW * _mV + _pW * _pV)
          + _nz * (_mW * _mW + _pW * _pW) + maxvel * (_mW - _pW));
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

  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    const DG_FP _nx = nx[1];
    const DG_FP _ny = ny[1];
    const DG_FP _nz = nz[1];
    const DG_FP _fscale = fscale[1];

    const DG_FP _mU = mU[i];
    const DG_FP _mV = mV[i];
    const DG_FP _mW = mW[i];
    const DG_FP _pU = pU[i];
    const DG_FP _pV = pV[i];
    const DG_FP _pW = pW[i];

    const DG_FP velL = _nx * _mU + _ny * _mV + _nz * _mW;
    const DG_FP velR = _nx * _pU + _ny * _pV + _nz * _pW;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[1][fIndCubR + i] = 0.5 * _fscale * (_nx * (_mU * _mU + _pU * _pU) + _ny * (_mU * _mV + _pU * _pV)
          + _nz * (_mU * _mW + _pU * _pW) + maxvel * (_mU - _pU));
    fV[1][fIndCubR + i] = 0.5 * _fscale * (_nx * (_mV * _mU + _pV * _pU) + _ny * (_mV * _mV + _pV * _pV)
          + _nz * (_mV * _mW + _pV * _pW) + maxvel * (_mV - _pV));
    fW[1][fIndCubR + i] = 0.5 * _fscale * (_nx * (_mW * _mU + _pW * _pU) + _ny * (_mW * _mV + _pW * _pV)
          + _nz * (_mW * _mW + _pW * _pW) + maxvel * (_mW - _pW));
  }
}
