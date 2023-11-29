inline void ls_advec_3d_oi_bc(const int *faceNum, const int *bc_type, const DG_FP *x,
                           const DG_FP *y, const DG_FP *z, const DG_FP *nx,
                           const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                           const DG_FP *val, const DG_FP *u, const DG_FP *v, const DG_FP *w,
                           DG_FP *flux) {
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;
  const int fIndCub = *faceNum * DG_CUB_SURF_3D_NP;

  DG_FP tmp_pVal[DG_NPF];
  if(*bc_type == BC_TYPE_NO_SLIP || *bc_type == BC_TYPE_SLIP || *bc_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      tmp_pVal[i] = val[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      tmp_pVal[i] = ps3d_custom_bc_get_ls(*bc_type, x[fmask_ind], y[fmask_ind], z[fmask_ind], val[fmask_ind]);
    }
  }

  DG_FP mU[DG_CUB_SURF_3D_NP], mV[DG_CUB_SURF_3D_NP], mW[DG_CUB_SURF_3D_NP];
  DG_FP mVal[DG_CUB_SURF_3D_NP], pVal[DG_CUB_SURF_3D_NP];
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    mU[i] = 0.0; mV[i] = 0.0; mW[i] = 0.0;
    mVal[i] = 0.0; pVal[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];

    const DG_FP _uL = u[fmask_ind];
    const DG_FP _vL = v[fmask_ind];
    const DG_FP _wL = w[fmask_ind];
    const DG_FP _valL = val[fmask_ind];
    const DG_FP _valR = tmp_pVal[i];

    for(int j = 0; j < DG_CUB_SURF_3D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCub + j, fInd + i, DG_CUB_SURF_3D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf3d_Interp_kernel[ind];
      mU[j] += mat_val * _uL;
      mV[j] += mat_val * _vL;
      mW[j] += mat_val * _wL;
      mVal[j] += mat_val * _valL;
      pVal[j] += mat_val * _valR;
    }
  }

  const DG_FP _nx = *nx;
  const DG_FP _ny = *ny;
  const DG_FP _nz = *nz;
  const DG_FP _fscale = *fscale;
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    const DG_FP _mU = mU[i];
    const DG_FP _mV = mV[i];
    const DG_FP _mW = mW[i];
    const DG_FP _mVal = mVal[i];
    const DG_FP _pVal = pVal[i];

    const DG_FP flux0 = _nx * _mU + _ny * _mV + _nz * _mW;
    const DG_FP flux1 = 0.5 * (_mVal + _pVal);
    const DG_FP flux2 = fabs(flux0);
    const DG_FP flux3 = _mVal - _pVal;

    flux[fIndCub + i] = _fscale * (flux0 * flux1 + 0.5 * flux2 * flux3);
  }
}
