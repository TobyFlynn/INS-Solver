inline void ins_2d_advec_sc_rhs_oi_2(const DG_FP *t, const int *bedge_type,
                        const int *bedgeNum, const DG_FP *nx, const DG_FP *ny,
                        const DG_FP *fscale, const DG_FP *x, const DG_FP *y,
                        const DG_FP *us, const DG_FP *vs, const DG_FP *ub,
                        const DG_FP *vb, DG_FP *fU, DG_FP *fV) {
  const int edge = *bedgeNum;
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];
  const int fInd = edge * DG_NPF;
  const int fIndCub = edge * DG_CUB_SURF_2D_NP;

  DG_FP usR[DG_NPF], vsR[DG_NPF], ubR[DG_NPF], vbR[DG_NPF];
  if(*bedge_type == BC_TYPE_NO_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      usR[i] = 0.0;
      vsR[i] = 0.0;
      ubR[i] = 0.0;
      vbR[i] = 0.0;
    }
  } else if (*bedge_type == BC_TYPE_SLIP) {
    DG_FP tangent_x = *ny;
    DG_FP tangent_y = -*nx;
    DG_FP tangent_mag = sqrt(tangent_x * tangent_x + tangent_y * tangent_y);
    tangent_x /= tangent_mag;
    tangent_y /= tangent_mag;
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      const DG_FP dot_s = us[fmask_ind] * tangent_x + vs[fmask_ind] * tangent_y;
      const DG_FP dot_b = ub[fmask_ind] * tangent_x + vb[fmask_ind] * tangent_y;
      usR[i] = dot_s * tangent_x;
      vsR[i] = dot_s * tangent_y;
      ubR[i] = dot_b * tangent_x;
      vbR[i] = dot_b * tangent_y;
    }
  } else if(*bedge_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      usR[i] = us[fmask_ind];
      vsR[i] = vs[fmask_ind];
      ubR[i] = ub[fmask_ind];
      vbR[i] = vb[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      ps2d_custom_bc_get_vel(*bedge_type, *t, x[fmask_ind], y[fmask_ind], *nx, *ny, us[fmask_ind], vs[fmask_ind], usR[i], vsR[i]);
      ps2d_custom_bc_get_vel(*bedge_type, *t, x[fmask_ind], y[fmask_ind], *nx, *ny, ub[fmask_ind], vb[fmask_ind], ubR[i], vbR[i]);
    }
  }

  DG_FP mUs[DG_CUB_SURF_2D_NP], mVs[DG_CUB_SURF_2D_NP];
  DG_FP mUb[DG_CUB_SURF_2D_NP], mVb[DG_CUB_SURF_2D_NP];
  DG_FP pUs[DG_CUB_SURF_2D_NP], pVs[DG_CUB_SURF_2D_NP];
  DG_FP pUb[DG_CUB_SURF_2D_NP], pVb[DG_CUB_SURF_2D_NP];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    mUs[i] = 0.0; mVs[i] = 0.0;
    mUb[i] = 0.0; mVb[i] = 0.0;
    pUs[i] = 0.0; pVs[i] = 0.0;
    pUb[i] = 0.0; pVb[i] = 0.0;
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];

    DG_FP _usL = us[fmask_ind];
    DG_FP _vsL = vs[fmask_ind];
    DG_FP _ubL = ub[fmask_ind];
    DG_FP _vbL = vb[fmask_ind];
    DG_FP _usR = usR[i];
    DG_FP _vsR = vsR[i];
    DG_FP _ubR = ubR[i];
    DG_FP _vbR = vbR[i];

    for(int j = 0; j < DG_CUB_SURF_2D_NP; j++) {
      const int ind = DG_MAT_IND(fIndCub + j, fInd + i, DG_CUB_SURF_2D_NP * DG_NUM_FACES, DG_NPF * DG_NUM_FACES);
      const DG_FP mat_val = dg_cubSurf2d_Interp_kernel[ind];
      mUs[j] += mat_val * _usL;
      mVs[j] += mat_val * _vsL;
      mUb[j] += mat_val * _ubL;
      mVb[j] += mat_val * _vbL;
      pUs[j] += mat_val * _usR;
      pVs[j] += mat_val * _vsR;
      pUb[j] += mat_val * _ubR;
      pVb[j] += mat_val * _vbR;
    }
  }

  const DG_FP _nx = *nx;
  const DG_FP _ny = *ny;
  const DG_FP _fscale = *fscale;
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    const DG_FP _mUs = mUs[i];
    const DG_FP _mVs = mVs[i];
    const DG_FP _mUb = mUb[i];
    const DG_FP _mVb = mVb[i];
    const DG_FP _pUs = pUs[i];
    const DG_FP _pVs = pVs[i];
    const DG_FP _pUb = pUb[i];
    const DG_FP _pVb = pVb[i];

    const DG_FP velL = _nx * _mUb + _ny * _mVb;
    const DG_FP velR = _nx * _pUb + _ny * _pVb;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    fU[fIndCub + i] += 0.5 * _fscale * (_nx * (_mUb * _mUs + _pUb * _pUs) + _ny * (_mUb * _mVs + _pUb * _pVs)
          + maxvel * (_mUs - _pUs));
    fV[fIndCub + i] += 0.5 * _fscale * (_nx * (_mVb * _mUs + _pVb * _pUs) + _ny * (_mVb * _mVs + _pVb * _pVs)
          + maxvel * (_mVs - _pVs));
  }
}
