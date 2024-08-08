inline void ins_2d_advec_oi_2(const DG_FP *t, const int *bedge_type,
                        const int *bedgeNum, const DG_FP *nx, const DG_FP *ny,
                        const DG_FP *fscale, const DG_FP *x, const DG_FP *y,
                        const DG_FP *u, const DG_FP *v, DG_FP *fU, DG_FP *fV) {
  const int edge = *bedgeNum;
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];
  const int fInd = edge * DG_NPF;
  const int fIndCub = edge * DG_CUB_SURF_2D_NP;

  DG_FP uR[DG_NPF], vR[DG_NPF];
  if(*bedge_type == BC_TYPE_NO_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      uR[i] = 0.0;
      vR[i] = 0.0;
    }
  } else if (*bedge_type == BC_TYPE_SLIP) {
    DG_FP tangent_x = *ny;
    DG_FP tangent_y = -*nx;
    DG_FP tangent_mag = sqrt(tangent_x * tangent_x + tangent_y * tangent_y);
    tangent_x /= tangent_mag;
    tangent_y /= tangent_mag;
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      const DG_FP dot = u[fmask_ind] * tangent_x + v[fmask_ind] * tangent_y;
      uR[i] = dot * tangent_x;
      vR[i] = dot * tangent_y;
    }
  } else if(*bedge_type == BC_TYPE_SLIP_X) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      uR[i] = 0.0;
      vR[i] = v[fmask_ind];
    }
  } else if(*bedge_type == BC_TYPE_SLIP_Y) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      uR[i] = u[fmask_ind];
      vR[i] = 0.0;
    }
  } else if(*bedge_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      uR[i] = u[fmask_ind];
      vR[i] = v[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      ps2d_custom_bc_get_vel(*bedge_type, *t, x[fmask_ind], y[fmask_ind], *nx,
                             *ny, u[fmask_ind], v[fmask_ind], uR[i], vR[i]);
    }
  }

  DG_FP mU[DG_CUB_SURF_2D_NP], mV[DG_CUB_SURF_2D_NP];
  DG_FP pU[DG_CUB_SURF_2D_NP], pV[DG_CUB_SURF_2D_NP];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    mU[i] = 0.0; mV[i] = 0.0;
    pU[i] = 0.0; pV[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];

    DG_FP _uL = u[fmask_ind];
    DG_FP _vL = v[fmask_ind];
    DG_FP _uR = uR[i];
    DG_FP _vR = vR[i];

    for(int j = 0; j < DG_CUB_SURF_2D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCub + j, fInd + i, DG_CUB_SURF_2D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf2d_Interp_kernel[ind];
      mU[j] += mat_val * _uL;
      mV[j] += mat_val * _vL;
      pU[j] += mat_val * _uR;
      pV[j] += mat_val * _vR;
    }
  }

  const DG_FP _nx = *nx;
  const DG_FP _ny = *ny;
  const DG_FP _fscale = *fscale;
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    const DG_FP _mU = mU[i];
    const DG_FP _mV = mV[i];
    const DG_FP _pU = pU[i];
    const DG_FP _pV = pV[i];

    const DG_FP velL = _nx * _mU + _ny * _mV;
    const DG_FP velR = _nx * _pU + _ny * _pV;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[fIndCub + i] += 0.5 * _fscale * (_nx * (_mU * _mU + _pU * _pU) + _ny * (_mU * _mV + _pU * _pV)
          + maxvel * (_mU - _pU));
    fV[fIndCub + i] += 0.5 * _fscale * (_nx * (_mV * _mU + _pV * _pU) + _ny * (_mV * _mV + _pV * _pV)
          + maxvel * (_mV - _pV));
  }
}
