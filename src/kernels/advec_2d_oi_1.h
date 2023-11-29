inline void advec_2d_oi_1(const int *faceNum, const bool *reverse,
                          const DG_FP *nx, const DG_FP *ny, const DG_FP *fscale,
                          const DG_FP **u, const DG_FP **v, const DG_FP **val,
                          DG_FP **flux) {
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
  DG_FP mVal[DG_CUB_SURF_2D_NP], pVal[DG_CUB_SURF_2D_NP];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    mU[i] = 0.0; mV[i] = 0.0;
    mVal[i] = 0.0; pVal[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];

    const DG_FP _uL = u[0][fmaskL_ind];
    const DG_FP _vL = v[0][fmaskL_ind];
    const DG_FP _valL = val[0][fmaskL_ind];
    const DG_FP _valR = val[1][fmaskR_ind];

    for(int j = 0; j < DG_CUB_SURF_2D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCubL + j, fIndL + i, DG_CUB_SURF_2D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf2d_Interp_kernel[ind];
      mU[j] += mat_val * _uL;
      mV[j] += mat_val * _vL;
      mVal[j] += mat_val * _valL;
      pVal[j] += mat_val * _valR;
    }
  }

  const DG_FP _nxL = nx[0];
  const DG_FP _nyL = ny[0];
  const DG_FP _fscaleL = fscale[0];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    const DG_FP _mU = mU[i];
    const DG_FP _mV = mV[i];
    const DG_FP _mVal = mVal[i];
    const DG_FP _pVal = pVal[i];

    const DG_FP flux0 = _nxL * _mU + _nyL * _mV;
    const DG_FP flux1 = 0.5 * (_mVal + _pVal);
    const DG_FP flux2 = fabs(flux0);
    const DG_FP flux3 = _mVal - _pVal;

    flux[0][fIndCubL + i] = _fscaleL * (flux0 * flux1 + 0.5 * flux2 * flux3);
  }

  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    mU[i] = 0.0; mV[i] = 0.0;
    mVal[i] = 0.0; pVal[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = rev ? fmaskL[DG_NPF - i - 1] : fmaskL[i];
    const int fmaskR_ind = fmaskR[i];

    const DG_FP _uR = u[1][fmaskR_ind];
    const DG_FP _vR = v[1][fmaskR_ind];
    const DG_FP _valL = val[0][fmaskL_ind];
    const DG_FP _valR = val[1][fmaskR_ind];

    for(int j = 0; j < DG_CUB_SURF_2D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCubR + j, fIndR + i, DG_CUB_SURF_2D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf2d_Interp_kernel[ind];
      mU[j] += mat_val * _uR;
      mV[j] += mat_val * _vR;
      mVal[j] += mat_val * _valR;
      pVal[j] += mat_val * _valL;
    }
  }

  const DG_FP _nxR = nx[1];
  const DG_FP _nyR = ny[1];
  const DG_FP _fscaleR = fscale[1];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    const DG_FP _mU = mU[i];
    const DG_FP _mV = mV[i];
    const DG_FP _mVal = mVal[i];
    const DG_FP _pVal = pVal[i];

    const DG_FP flux0 = _nxR * _mU + _nyR * _mV;
    const DG_FP flux1 = 0.5 * (_mVal + _pVal);
    const DG_FP flux2 = fabs(flux0);
    const DG_FP flux3 = _mVal - _pVal;

    flux[1][fIndCubR + i] = _fscaleR * (flux0 * flux1 + 0.5 * flux2 * flux3);
  }
}
