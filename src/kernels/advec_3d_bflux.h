inline void advec_3d_bflux(const int *faceNum, const int *bc_type, const DG_FP *x,
                           const DG_FP *y, const DG_FP *z, const DG_FP *nx,
                           const DG_FP *ny, const DG_FP *nz,
                           const DG_FP *fscale, const DG_FP *val,
                           const DG_FP *u, const DG_FP *v, const DG_FP *w,
                           DG_FP *flux) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * DG_NPF];

  if(*bc_type == LW_INFLOW_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      int find = *faceNum * DG_NPF + i;
      DG_FP flux0 = *nx * u[fmask[i]] + *ny * v[fmask[i]] + *nz * w[fmask[i]];
      DG_FP flux1 = sqrt(x[fmask[i]] * x[fmask[i]] + y[fmask[i]] * y[fmask[i]] + z[fmask[i]] * z[fmask[i]]) - LW_INLET_RADIUS;

      flux[find] += *fscale * flux0 * flux1;
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      int find = *faceNum * DG_NPF + i;
      DG_FP flux0 = *nx * u[fmask[i]] + *ny * v[fmask[i]] + *nz * w[fmask[i]];
      DG_FP flux1 = val[fmask[i]];

      flux[find] += *fscale * flux0 * flux1;
    }
  }
}
