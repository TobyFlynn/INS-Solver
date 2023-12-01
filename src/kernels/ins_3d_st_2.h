inline void ins_3d_st_2(const DG_FP *alpha_, const int *faceNum, const int *fmaskL_corrected,
                        const int *fmaskR_corrected, const DG_FP *nx,
                        const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                        const DG_FP **s, DG_FP **fX, DG_FP **fY, DG_FP **fZ) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];
  const int fIndL = faceNum[0] * DG_NPF;
  const int fIndR = faceNum[1] * DG_NPF;
  const int fIndCubL = faceNum[0] * DG_CUB_SURF_3D_NP;
  const int fIndCubR = faceNum[1] * DG_CUB_SURF_3D_NP;
  const DG_FP alpha = *alpha_;
  const DG_FP PI = 3.141592653589793238463;

  DG_FP mS[DG_CUB_SURF_3D_NP];
  DG_FP pS[DG_CUB_SURF_3D_NP];
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    mS[i] = 0.0;
    pS[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];

    const DG_FP _sL = s[0][fmaskL_ind];
    const DG_FP _sR = s[1][fmaskR_ind];

    for(int j = 0; j < DG_CUB_SURF_3D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCubL + j, fIndL + i, DG_CUB_SURF_3D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf3d_Interp_kernel[ind];
      mS[j] += mat_val * _sL;
      pS[j] += mat_val * _sR;
    }
  }

  const DG_FP _nxL = nx[0];
  const DG_FP _nyL = ny[0];
  const DG_FP _nzL = nz[0];
  const DG_FP _fscaleL = fscale[0];
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    const DG_FP _mS = 0.5 * tanh(PI * mS[i] / alpha);
    const DG_FP _pS = 0.5 * tanh(PI * pS[i] / alpha);
    const DG_FP avg = 0.5 * (_mS + _pS);

    fX[0][fIndCubL + i] = _fscaleL * _nxL * avg;
    fY[0][fIndCubL + i] = _fscaleL * _nyL * avg;
    fZ[0][fIndCubL + i] = _fscaleL * _nzL * avg;
  }

  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    mS[i] = 0.0;
    pS[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL_corrected[i];
    const int fmaskR_ind = fmaskR[i];

    const DG_FP _sL = s[0][fmaskL_ind];
    const DG_FP _sR = s[1][fmaskR_ind];

    for(int j = 0; j < DG_CUB_SURF_3D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCubR + j, fIndR + i, DG_CUB_SURF_3D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf3d_Interp_kernel[ind];
      mS[j] += mat_val * _sR;
      pS[j] += mat_val * _sL;
    }
  }

  const DG_FP _nxR = nx[1];
  const DG_FP _nyR = ny[1];
  const DG_FP _nzR = nz[1];
  const DG_FP _fscaleR = fscale[1];
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    const DG_FP _mS = 0.5 * tanh(PI * mS[i] / alpha);
    const DG_FP _pS = 0.5 * tanh(PI * pS[i] / alpha);
    const DG_FP avg = 0.5 * (_mS + _pS);

    fX[1][fIndCubR + i] = _fscaleR * _nxR * avg;
    fY[1][fIndCubR + i] = _fscaleR * _nyR * avg;
    fZ[1][fIndCubR + i] = _fscaleR * _nzR * avg;
  }
}
