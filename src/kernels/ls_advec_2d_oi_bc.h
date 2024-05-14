inline void ls_advec_2d_oi_bc(const int *bedge_type, const int *bedgeNum,
                          const DG_FP *nx, const DG_FP *ny, const DG_FP *fscale,
                          const DG_FP *x, const DG_FP *y, const DG_FP *u,
                          const DG_FP *v, const DG_FP *val, DG_FP *flux) {
  const int edge = bedgeNum[0];
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];
  const int fInd = edge * DG_NPF;
  const int fIndCub = edge * DG_CUB_SURF_2D_NP;

  DG_FP tmp_valP[DG_NPF];
  if(*bedge_type == BC_TYPE_NO_SLIP || *bedge_type == BC_TYPE_SLIP || *bedge_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      tmp_valP[i] = val[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      tmp_valP[i] = ps2d_custom_bc_get_ls(*bedge_type, x[fmask_ind], y[fmask_ind], val[fmask_ind]);
    }
  }

  DG_FP mU[DG_CUB_SURF_2D_NP], mV[DG_CUB_SURF_2D_NP];
  DG_FP mVal[DG_CUB_SURF_2D_NP], pVal[DG_CUB_SURF_2D_NP];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    mU[i] = 0.0; mV[i] = 0.0;
    mVal[i] = 0.0; pVal[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];

    const DG_FP _uL = u[fmask_ind];
    const DG_FP _vL = v[fmask_ind];
    const DG_FP _valL = val[fmask_ind];
    const DG_FP _valR = tmp_valP[i];

    for(int j = 0; j < DG_CUB_SURF_2D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCub + j, fInd + i, DG_CUB_SURF_2D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf2d_Interp_kernel[ind];
      mU[j] += mat_val * _uL;
      mV[j] += mat_val * _vL;
      mVal[j] += mat_val * _valL;
      pVal[j] += mat_val * _valR;
    }
  }

  const DG_FP _nx = *nx;
  const DG_FP _ny = *ny;
  const DG_FP _fscale = *fscale;
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    const DG_FP _mU = mU[i];
    const DG_FP _mV = mV[i];
    const DG_FP _mVal = mVal[i];
    const DG_FP _pVal = pVal[i];

    const DG_FP flux0 = _nx * _mU + _ny * _mV;
    const DG_FP flux1 = 0.5 * (_mVal + _pVal);
    const DG_FP flux2 = fabs(flux0);
    const DG_FP flux3 = _mVal - _pVal;

    flux[fIndCub + i] += _fscale * (flux0 * flux1 + 0.5 * flux2 * flux3);
  }
}
