inline void ins_2d_proj_7(const DG_FP *g0, const DG_FP *t, const int *bedge_type,
                          const int *bedgeNum, const DG_FP *nx, const DG_FP *ny,
                          const DG_FP *sJ, const DG_FP *x, const DG_FP *y,
                          const DG_FP *pen, const DG_FP *u, const DG_FP *v,
                          DG_FP *f0, DG_FP *f1) {
  const int edgeNum = bedgeNum[0];
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeNum * DG_NPF];
  const int fInd = edgeNum * DG_NPF;

  // Set boundary velocities
  DG_FP pU[DG_NPF], pV[DG_NPF];
  if(*bedge_type == BC_TYPE_NO_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      pU[i] = 0.0;
      pV[i] = 0.0;
    }
  } else if(*bedge_type == BC_TYPE_SLIP_X) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      pU[i] = 0.0;
      pV[i] = v[fmask_ind];
    }
  } else if(*bedge_type == BC_TYPE_SLIP_Y) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      pU[i] = u[fmask_ind];
      pV[i] = 0.0;
    }
  } else if(*bedge_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      pU[i] = u[fmask_ind];
      pV[i] = v[fmask_ind];
    }
  } else {
    const DG_FP _g0 = *g0;
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      ps2d_custom_bc_get_vel(*bedge_type, *t, x[fmask_ind], y[fmask_ind], nx[0], ny[0], u[fmask_ind] / _g0, v[fmask_ind] / _g0, pU[i], pV[i]);
      pU[i] *= _g0;
      pV[i] *= _g0;
    }
  }

  const DG_FP pen_f = *pen;
  const DG_FP _sJ = *sJ;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    f0[fInd + i] += _sJ * pen_f * (u[fmask_ind] - pU[i]);
    f1[fInd + i] += _sJ * pen_f * (v[fmask_ind] - pV[i]);
  }
}
