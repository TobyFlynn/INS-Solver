inline void ins_3d_advec_sc_rhs_oi_2(const DG_FP *t, const int *bc_type,
                          const int *faceNum, const DG_FP *nx, const DG_FP *ny,
                          const DG_FP *nz, const DG_FP *fscale, const DG_FP *x,
                          const DG_FP *y, const DG_FP *z, const DG_FP *us,
                          const DG_FP *vs, const DG_FP *ws, const DG_FP *ub,
                          const DG_FP *vb, const DG_FP *wb, DG_FP *fU,
                          DG_FP *fV, DG_FP *fW) {
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;
  const int fIndCub = *faceNum * DG_CUB_SURF_3D_NP;

  DG_FP mUs[DG_CUB_SURF_3D_NP], mVs[DG_CUB_SURF_3D_NP], mWs[DG_CUB_SURF_3D_NP];
  DG_FP mUb[DG_CUB_SURF_3D_NP], mVb[DG_CUB_SURF_3D_NP], mWb[DG_CUB_SURF_3D_NP];
  DG_FP pUs[DG_CUB_SURF_3D_NP], pVs[DG_CUB_SURF_3D_NP], pWs[DG_CUB_SURF_3D_NP];
  DG_FP pUb[DG_CUB_SURF_3D_NP], pVb[DG_CUB_SURF_3D_NP], pWb[DG_CUB_SURF_3D_NP];
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    mUs[i] = 0.0; mVs[i] = 0.0; mWs[i] = 0.0;
    mUb[i] = 0.0; mVb[i] = 0.0; mWb[i] = 0.0;
    pUs[i] = 0.0; pVs[i] = 0.0; pWs[i] = 0.0;
    pUb[i] = 0.0; pVb[i] = 0.0; pWb[i] = 0.0;
  }

  DG_FP bcUs[DG_NPF], bcVs[DG_NPF], bcWs[DG_NPF];
  DG_FP bcUb[DG_NPF], bcVb[DG_NPF], bcWb[DG_NPF];
  if(*bc_type == BC_TYPE_NO_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      bcUs[i] = 0.0;
      bcVs[i] = 0.0;
      bcWs[i] = 0.0;
      bcUb[i] = 0.0;
      bcVb[i] = 0.0;
      bcWb[i] = 0.0;
    }
  } else if(*bc_type == BC_TYPE_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      // S
      const DG_FP mags = sqrt(us[fmask_ind] * us[fmask_ind] + vs[fmask_ind] * vs[fmask_ind] + ws[fmask_ind] * ws[fmask_ind]);
      const DG_FP dots = *nx * us[fmask_ind] + *ny * vs[fmask_ind] + *nz * ws[fmask_ind];
      bcUs[i] = us[fmask_ind] - dots * *nx;
      bcVs[i] = vs[fmask_ind] - dots * *ny;
      bcWs[i] = ws[fmask_ind] - dots * *nz;
      const DG_FP mag2s = sqrt(bcUs[i] * bcUs[i] + bcVs[i] * bcVs[i] + bcWs[i] * bcWs[i]);
      const DG_FP mag_factors = fabs(mags) < 1e-8 || fabs(mag2s) < 1e-8 ? 1.0 : mags / mag2s;
      bcUs[i] *= mag_factors;
      bcVs[i] *= mag_factors;
      bcWs[i] *= mag_factors;
      // B
      const DG_FP magb = sqrt(ub[fmask_ind] * ub[fmask_ind] + vb[fmask_ind] * vb[fmask_ind] + wb[fmask_ind] * wb[fmask_ind]);
      const DG_FP dotb = *nx * ub[fmask_ind] + *ny * vb[fmask_ind] + *nz * wb[fmask_ind];
      bcUb[i] = ub[fmask_ind] - dotb * *nx;
      bcVs[i] = vb[fmask_ind] - dotb * *ny;
      bcWb[i] = wb[fmask_ind] - dotb * *nz;
      const DG_FP mag2b = sqrt(bcUb[i] * bcUb[i] + bcVs[i] * bcVs[i] + bcWb[i] * bcWb[i]);
      const DG_FP mag_factorb = fabs(magb) < 1e-8 || fabs(mag2b) < 1e-8 ? 1.0 : magb / mag2b;
      bcUb[i] *= mag_factorb;
      bcVs[i] *= mag_factorb;
      bcWb[i] *= mag_factorb;
    }
  } else if(*bc_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      bcUs[i] = us[fmask_ind];
      bcVs[i] = vs[fmask_ind];
      bcWs[i] = ws[fmask_ind];
      bcUb[i] = ub[fmask_ind];
      bcVb[i] = vb[fmask_ind];
      bcWb[i] = wb[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      ps3d_custom_bc_get_vel(*bc_type, *t, x[fmask_ind], y[fmask_ind], z[fmask_ind],
                             *nx, *ny, *nz, us[fmask_ind], vs[fmask_ind], ws[fmask_ind],
                             bcUs[i], bcVs[i], bcWs[i]);
      ps3d_custom_bc_get_vel(*bc_type, *t, x[fmask_ind], y[fmask_ind], z[fmask_ind],
                             *nx, *ny, *nz, ub[fmask_ind], vb[fmask_ind], wb[fmask_ind],
                             bcUb[i], bcVb[i], bcWb[i]);
    }
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];

    DG_FP _usL = us[fmask_ind];
    DG_FP _vsL = vs[fmask_ind];
    DG_FP _wsL = ws[fmask_ind];
    DG_FP _ubL = ub[fmask_ind];
    DG_FP _vbL = vb[fmask_ind];
    DG_FP _wbL = wb[fmask_ind];
    DG_FP _usR = bcUs[i];
    DG_FP _vsR = bcVs[i];
    DG_FP _wsR = bcWs[i];
    DG_FP _ubR = bcUb[i];
    DG_FP _vbR = bcVb[i];
    DG_FP _wbR = bcWb[i];

    for(int j = 0; j < DG_CUB_SURF_3D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCub + j, fInd + i, DG_CUB_SURF_3D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf3d_Interp_kernel[ind];
      mUs[j] += mat_val * _usL;
      mVs[j] += mat_val * _vsL;
      mWs[j] += mat_val * _wsL;
      mUb[j] += mat_val * _ubL;
      mVb[j] += mat_val * _vbL;
      mWb[j] += mat_val * _wbL;
      pUs[j] += mat_val * _usR;
      pVs[j] += mat_val * _vsR;
      pWs[j] += mat_val * _wsR;
      pUb[j] += mat_val * _ubR;
      pVb[j] += mat_val * _vbR;
      pWb[j] += mat_val * _wbR;
    }
  }

  const DG_FP _nx = *nx;
  const DG_FP _ny = *ny;
  const DG_FP _nz = *nz;
  const DG_FP _fscale = *fscale;
  for(int i = 0; i < DG_CUB_SURF_3D_NP; i++) {
    const DG_FP _mUs = mUs[i];
    const DG_FP _mVs = mVs[i];
    const DG_FP _mWs = mWs[i];
    const DG_FP _pUs = pUs[i];
    const DG_FP _pVs = pVs[i];
    const DG_FP _pWs = pWs[i];
    const DG_FP _mUb = mUb[i];
    const DG_FP _mVb = mVb[i];
    const DG_FP _mWb = mWb[i];
    const DG_FP _pUb = pUb[i];
    const DG_FP _pVb = pVb[i];
    const DG_FP _pWb = pWb[i];

    // Check whether it shoudl be s instead of b here (or the max of them)
    const DG_FP velL = _nx * _mUb + _ny * _mVb + _nz * _mWb;
    const DG_FP velR = _nx * _pUb + _ny * _pVb + _nz * _pWb;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[fIndCub + i] += 0.5 * _fscale * (_nx * (_mUb * _mUs + _pUb * _pUs) + _ny * (_mUb * _mVs + _pUb * _pVs)
          + _nz * (_mUb * _mWs + _pUb * _pWs) + maxvel * (_mUs - _pUs));
    fV[fIndCub + i] += 0.5 * _fscale * (_nx * (_mVb * _mUs + _pVb * _pUs) + _ny * (_mVb * _mVs + _pVb * _pVs)
          + _nz * (_mVb * _mWs + _pVb * _pWs) + maxvel * (_mVs - _pVs));
    fW[fIndCub + i] += 0.5 * _fscale * (_nx * (_mWb * _mUs + _pWb * _pUs) + _ny * (_mWb * _mVs + _pWb * _pVs)
          + _nz * (_mWb * _mWs + _pWb * _pWs) + maxvel * (_mWs - _pWs));
  }
}
