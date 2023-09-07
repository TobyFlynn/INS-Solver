inline void advec_2d_bflux(const int *edgeNum, const DG_FP *nx,
                           const DG_FP *ny, const DG_FP *fscale,
                           const DG_FP *val, const DG_FP *u, const DG_FP *v,
                           DG_FP *flux) {
  // Work out which edge for each element
  const int edge = edgeNum[0];
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];

  const DG_FP _nx = nx[0];
  const DG_FP _ny = ny[0];
  const DG_FP _fscale = fscale[0];
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];

    DG_FP flux0 = _nx * u[fmask_ind] + _ny * v[fmask_ind];
    DG_FP flux1 = val[fmask_ind];
    DG_FP flux2 = fabs(flux0);
    DG_FP flux3 = 0.0;

    const int outL_ind = edge * DG_NPF + i;

    flux[outL_ind] += _fscale * (flux0 * flux1 + 0.5 * flux2 * flux3);
  }
}
