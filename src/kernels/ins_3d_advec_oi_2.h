inline void ins_3d_advec_oi_2(const DG_FP *t, const int *bc_type,
                          const int *faceNum, const DG_FP *nx, const DG_FP *ny,
                          const DG_FP *nz, const DG_FP *fscale, const DG_FP *x,
                          const DG_FP *y, const DG_FP *z, const DG_FP *u,
                          const DG_FP *v, const DG_FP *w, DG_FP *fU, DG_FP *fV,
                          DG_FP *fW) {
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;
  const int fIndCub = *faceNum * DG_CUB_SURF_3D_NP;

  DG_FP uR[DG_NPF], vR[DG_NPF], wR[DG_NPF];
  if(*bc_type == BC_TYPE_NO_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      uR[i] = 0.0;
      vR[i] = 0.0;
      wR[i] = 0.0;
    }
  } else if(*bc_type == BC_TYPE_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      const DG_FP mag = sqrt(u[fmask_ind] * u[fmask_ind] + v[fmask_ind] * v[fmask_ind] + w[fmask_ind] * w[fmask_ind]);
      const DG_FP dot = *nx * u[fmask_ind] + *ny * v[fmask_ind] + *nz * w[fmask_ind];
      uR[i] = u[fmask_ind] - dot * *nx;
      vR[i] = v[fmask_ind] - dot * *ny;
      wR[i] = w[fmask_ind] - dot * *nz;
      const DG_FP mag2 = sqrt(uR[fInd + i] * uR[fInd + i] + vR[fInd + i] * vR[fInd + i] + wR[fInd + i] * wR[fInd + i]);
      const DG_FP mag_factor = fabs(mag) < 1e-8 || fabs(mag2) < 1e-8 ? 1.0 : mag / mag2;
      uR[i] *= mag_factor;
      vR[i] *= mag_factor;
      wR[i] *= mag_factor;
    }
  } else if(*bc_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      uR[i] = u[fmask_ind];
      vR[i] = v[fmask_ind];
      wR[i] = w[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      ps3d_custom_bc_get_vel(*bc_type, *t, x[fmask_ind], y[fmask_ind], z[fmask_ind],
                             *nx, *ny, *nz, u[fmask_ind], v[fmask_ind], w[fmask_ind],
                             uR[i], vR[i], wR[i]);
    }
  }

  DG_FP mU[DG_CUB_SURF_3D_NP], mV[DG_CUB_SURF_3D_NP], mW[DG_CUB_SURF_3D_NP];
  DG_FP pU[DG_CUB_SURF_3D_NP], pV[DG_CUB_SURF_3D_NP], pW[DG_CUB_SURF_3D_NP];
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    mU[i] = 0.0; mV[i] = 0.0; mW[i] = 0.0;
    pU[i] = 0.0; pV[i] = 0.0; pW[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];

    DG_FP _uL = u[fmask_ind];
    DG_FP _vL = v[fmask_ind];
    DG_FP _wL = w[fmask_ind];
    DG_FP _uR = uR[i];
    DG_FP _vR = vR[i];
    DG_FP _wR = wR[i];

    for(int j = 0; j < DG_CUB_SURF_3D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCub + j, fInd + i, DG_CUB_SURF_3D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf3d_Interp_kernel[ind];
      mU[j] += mat_val * _uL;
      mV[j] += mat_val * _vL;
      mW[j] += mat_val * _wL;
      pU[j] += mat_val * _uR;
      pV[j] += mat_val * _vR;
      pW[j] += mat_val * _wR;
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
    const DG_FP _pU = pU[i];
    const DG_FP _pV = pV[i];
    const DG_FP _pW = pW[i];

    const DG_FP velL = _nx * _mU + _ny * _mV + _nz * _mW;
    const DG_FP velR = _nx * _pU + _ny * _pV + _nz * _pW;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[fIndCub + i] += 0.5 * _fscale * (_nx * (_mU * _mU + _pU * _pU) + _ny * (_mU * _mV + _pU * _pV)
          + _nz * (_mU * _mW + _pU * _pW) + maxvel * (_mU - _pU));
    fV[fIndCub + i] += 0.5 * _fscale * (_nx * (_mV * _mU + _pV * _pU) + _ny * (_mV * _mV + _pV * _pV)
          + _nz * (_mV * _mW + _pV * _pW) + maxvel * (_mV - _pV));
    fW[fIndCub + i] += 0.5 * _fscale * (_nx * (_mW * _mU + _pW * _pU) + _ny * (_mW * _mV + _pW * _pV)
          + _nz * (_mW * _mW + _pW * _pW) + maxvel * (_mW - _pW));
  }
}
