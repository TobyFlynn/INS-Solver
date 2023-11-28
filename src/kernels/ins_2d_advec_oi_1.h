inline void ins_2d_advec_oi_1(const int *faceNum, const bool *reverse,
                    const DG_FP *nx, const DG_FP *ny, const DG_FP *fscale,
                    const DG_FP **u, const DG_FP **v, DG_FP **fU, DG_FP **fV) {
  const bool rev = *reverse;
  const int edgeL = faceNum[0];
  const int edgeR = faceNum[1];
  const int *fmaskL = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeL * DG_NPF];
  const int *fmaskR = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeR * DG_NPF];
  const int fIndL = edgeL * DG_NPF;
  const int fIndR = edgeR * DG_NPF;
  const int fIndCubL = edgeL * DG_CUB_SURF_2D_NP;
  const int fIndCubR = edgeR * DG_CUB_SURF_2D_NP;

  DG_FP mU[DG_CUB_SURF_2D_NP], mV[DG_CUB_SURF_2D_NP];
  DG_FP pU[DG_CUB_SURF_2D_NP], pV[DG_CUB_SURF_2D_NP];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    mU[i] = 0.0; mV[i] = 0.0;
    pU[i] = 0.0; pV[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];

    DG_FP _uL = u[0][fmaskL_ind];
    DG_FP _vL = v[0][fmaskL_ind];
    DG_FP _uR = u[1][fmaskR_ind];
    DG_FP _vR = v[1][fmaskR_ind];

    for(int j = 0; j < DG_CUB_SURF_2D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCubL + j, fIndL + i, DG_CUB_SURF_2D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf2d_Interp_kernel[ind];
      mU[j] += mat_val * _uL;
      mV[j] += mat_val * _vL;
      pU[j] += mat_val * _uR;
      pV[j] += mat_val * _vR;
    }
  }

  const DG_FP _nxL = nx[0];
  const DG_FP _nyL = ny[0];
  const DG_FP _fscaleL = fscale[0];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    const DG_FP _mU = mU[i];
    const DG_FP _mV = mV[i];
    const DG_FP _pU = pU[i];
    const DG_FP _pV = pV[i];

    const DG_FP velL = _nxL * _mU + _nyL * _mV;
    const DG_FP velR = _nxL * _pU + _nyL * _pV;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[0][fIndCubL + i] = 0.5 * _fscaleL * (_nxL * (_mU * _mU + _pU * _pU) + _nyL * (_mU * _mV + _pU * _pV)
          + maxvel * (_mU - _pU));
    fV[0][fIndCubL + i] = 0.5 * _fscaleL * (_nxL * (_mV * _mU + _pV * _pU) + _nyL * (_mV * _mV + _pV * _pV)
          + maxvel * (_mV - _pV));
  }

  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    mU[i] = 0.0; mV[i] = 0.0;
    pU[i] = 0.0; pV[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskR_ind = fmaskR[i];
    const int fmaskL_ind = rev ? fmaskL[DG_NPF - i - 1] : fmaskL[i];

    DG_FP _uL = u[0][fmaskL_ind];
    DG_FP _vL = v[0][fmaskL_ind];
    DG_FP _uR = u[1][fmaskR_ind];
    DG_FP _vR = v[1][fmaskR_ind];

    for(int j = 0; j < DG_CUB_SURF_2D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCubR + j, fIndR + i, DG_CUB_SURF_2D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf2d_Interp_kernel[ind];
      pU[j] += mat_val * _uL;
      pV[j] += mat_val * _vL;
      mU[j] += mat_val * _uR;
      mV[j] += mat_val * _vR;
    }
  }

  const DG_FP _nxR = nx[1];
  const DG_FP _nyR = ny[1];
  const DG_FP _fscaleR = fscale[1];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    const DG_FP _mU = mU[i];
    const DG_FP _mV = mV[i];
    const DG_FP _pU = pU[i];
    const DG_FP _pV = pV[i];

    const DG_FP velL = _nxR * _mU + _nyR * _mV;
    const DG_FP velR = _nxR * _pU + _nyR * _pV;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[1][fIndCubR + i] = 0.5 * _fscaleR * (_nxR * (_mU * _mU + _pU * _pU) + _nyR * (_mU * _mV + _pU * _pV)
          + maxvel * (_mU - _pU));
    fV[1][fIndCubR + i] = 0.5 * _fscaleR * (_nxR * (_mV * _mU + _pV * _pU) + _nyR * (_mV * _mV + _pV * _pV)
          + maxvel * (_mV - _pV));
  }
}
