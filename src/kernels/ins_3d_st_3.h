inline void ins_3d_st_3(const DG_FP *alpha_, const int *bc_type, const int *faceNum,
                        const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                        const DG_FP *fscale, const DG_FP *s, DG_FP *fX,
                        DG_FP *fY, DG_FP *fZ) {
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;
  const int fIndCub = *faceNum * DG_CUB_SURF_3D_NP;
  const DG_FP alpha = *alpha_;
  const DG_FP PI = 3.141592653589793238463;

  // TODO treat boundary conditions differently
  DG_FP sR[DG_NPF];
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    sR[i] = s[fmask_ind];
  }

  DG_FP mS[DG_CUB_SURF_3D_NP];
  DG_FP pS[DG_CUB_SURF_3D_NP];
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    mS[i] = 0.0;
    pS[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    const DG_FP _sL = s[fmask_ind];
    const DG_FP _sR = sR[i];

    for(int j = 0; j < DG_CUB_SURF_3D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCub + j, fInd + i, DG_CUB_SURF_3D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf3d_Interp_kernel[ind];
      mS[j] += mat_val * _sL;
      pS[j] += mat_val * _sR;
    }
  }

  const DG_FP _nx = *nx;
  const DG_FP _ny = *ny;
  const DG_FP _nz = *nz;
  const DG_FP _fscale = *fscale;
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    const DG_FP _mS = 0.5 * tanh(PI * mS[i] / alpha);
    const DG_FP _pS = 0.5 * tanh(PI * pS[i] / alpha);
    const DG_FP avg = 0.5 * (_mS + _pS);

    fX[fIndCub + i] += _fscale * _nx * avg;
    fY[fIndCub + i] += _fscale * _ny * avg;
    fZ[fIndCub + i] += _fscale * _nz * avg;
  }
}
