inline void ls_advec_3d_bc(const int *faceNum, const int *bc_type, const DG_FP *x,
                           const DG_FP *y, const DG_FP *z, const DG_FP *nx,
                           const DG_FP *ny, const DG_FP *nz,
                           const DG_FP *fscale, const DG_FP *val,
                           const DG_FP *u, const DG_FP *v, const DG_FP *w,
                           DG_FP *flux) {
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];

  DG_FP sP[DG_NPF];
  if(*bc_type == BC_TYPE_NO_SLIP || *bc_type == BC_TYPE_SLIP || *bc_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      sP[i] = val[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      sP[i] = ps3d_custom_bc_get_ls(*bc_type, x[fmask_ind], y[fmask_ind], z[fmask_ind], val[fmask_ind]);
    }
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int find = *faceNum * DG_NPF + i;
    const int fmask_ind = fmask[i];
    DG_FP flux0 = *nx * u[fmask_ind] + *ny * v[fmask_ind] + *nz * w[fmask_ind];
    DG_FP flux1 = 0.5 * (val[fmask_ind] + sP[i]);
    DG_FP flux2 = fabs(flux0);
    DG_FP flux3 = val[fmask_ind] - sP[i];
    DG_FP flux4 = flux0 * flux1 + 0.5 * flux2 * flux3;

    flux[find] += *fscale * flux4;
  }
}
