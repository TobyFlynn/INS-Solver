inline void ls_advec_2d_bc(const int *bedge_type, const int *edgeNum, const DG_FP *nx, 
                           const DG_FP *ny, const DG_FP *fscale, const DG_FP *x, 
                           const DG_FP *y, const DG_FP *val, const DG_FP *u, 
                           const DG_FP *v, DG_FP *flux) {
  // Work out which edge for each element
  const int edge = edgeNum[0];
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];

  // Get BC value of LS
  DG_FP sP[DG_NPF];
  if(*bedge_type == BC_TYPE_NO_SLIP || *bedge_type == BC_TYPE_SLIP || *bedge_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      sP[i] = val[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      sP[i] = ps2d_custom_bc_get_ls(*bedge_type, x[fmask_ind], y[fmask_ind], val[fmask_ind]);
    }
  }

  const DG_FP _nx = nx[0];
  const DG_FP _ny = ny[0];
  const DG_FP _fscale = fscale[0];
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];

    DG_FP flux0 = _nx * u[fmask_ind] + _ny * v[fmask_ind];
    DG_FP flux1 = 0.5 * (val[fmask_ind] + sP[i]);
    DG_FP flux2 = fabs(flux0);
    DG_FP flux3 = val[fmask_ind] - sP[i];

    const int outL_ind = edge * DG_NPF + i;

    flux[outL_ind] += _fscale * (flux0 * flux1 + 0.5 * flux2 * flux3);
  }
}
