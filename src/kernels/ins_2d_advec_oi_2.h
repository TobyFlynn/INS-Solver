inline void ins_2d_advec_oi_2(const DG_FP *t, const int *bedge_type, const int *bedgeNum,
                              const DG_FP *nx, const DG_FP *ny, const DG_FP *x, const DG_FP *y,
                              const DG_FP *u, const DG_FP *v, DG_FP *uM, DG_FP *vM, DG_FP *uP,
                              DG_FP *vP) {
  const int edge = bedgeNum[0];
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];

  for(int i = 0; i < DG_NPF; i++) {
    const int find = fmask[i];
    uM[edge * DG_NPF + i] = u[find];
    vM[edge * DG_NPF + i] = v[find];
  }

  // Set boundary velocities
  if(*bedge_type == BC_TYPE_NO_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      uP[edge * DG_NPF + i] = 0.0;
      vP[edge * DG_NPF + i] = 0.0;
    }
  } else if(*bedge_type == BC_TYPE_SLIP_X) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      uP[edge * DG_NPF + i] = 0.0;
      vP[edge * DG_NPF + i] = v[fmask_ind];
    }
  } else if (*bedge_type == BC_TYPE_SLIP_Y) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      uP[edge * DG_NPF + i] = u[fmask_ind];
      vP[edge * DG_NPF + i] = 0.0;
    }
  } else if(*bedge_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      uP[edge * DG_NPF + i] = u[fmask_ind];
      vP[edge * DG_NPF + i] = v[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      ps2d_custom_bc_get_vel(*bedge_type, *t, x[fmask_ind], y[fmask_ind], nx[0], ny[0], u[fmask_ind], v[fmask_ind], uP[edge * DG_NPF + i], vP[edge * DG_NPF + i]);
    }
  }
}
