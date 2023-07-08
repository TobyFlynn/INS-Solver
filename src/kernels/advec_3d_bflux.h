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
      const int find = *faceNum * DG_NPF + i;
      const int fmask_ind = fmaskB[i];
      DG_FP flux0 = *nx * u[fmask_ind] + *ny * v[fmask_ind] + *nz * w[fmask_ind];
      DG_FP flux1 = 0.5 * (val[fmask_ind] + (-1.0));
      DG_FP flux2 = fabs(flux0);
      DG_FP flux3 = val[fmask_ind] - (-1.0);
      DG_FP flux4 = flux0 * flux1 + 0.5 * flux2 * flux3;

      flux[find] += *fscale * flux4;
      // flux[find] += *fscale * flux0 * (-1.0);
      // flux[find] += *fscale * flux0 * (val[fmaskB[i]] + 1.0);
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int find = *faceNum * DG_NPF + i;
      const int fmask_ind = fmaskB[i];
      DG_FP flux0 = *nx * u[fmask_ind] + *ny * v[fmask_ind] + *nz * w[fmask_ind];

      flux[find] += *fscale * flux0 * val[fmask_ind];
    }
  }
}
