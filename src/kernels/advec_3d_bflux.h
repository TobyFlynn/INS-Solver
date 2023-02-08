inline void advec_3d_bflux(const int *faceNum, const int *bc_type, const double *x,
                           const double *y, const double *z, const double *nx,
                           const double *ny, const double *nz,
                           const double *fscale, const double *val,
                           const double *u, const double *v, const double *w,
                           double *flux) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * DG_NPF];

  if(*bc_type == 0 || *bc_type == 5) {
    for(int i = 0; i < DG_NPF; i++) {
      int find = *faceNum * DG_NPF + i;
      double flux0 = *nx * u[fmask[i]] + *ny * v[fmask[i]] + *nz * w[fmask[i]];
      double flux1 = sqrt(y[fmask[i]] * y[fmask[i]] + z[fmask[i]] * z[fmask[i]]) - 0.05;

      flux[find] += *fscale * flux0 * flux1;
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      int find = *faceNum * DG_NPF + i;
      double flux0 = *nx * u[fmask[i]] + *ny * v[fmask[i]] + *nz * w[fmask[i]];
      double flux1 = val[fmask[i]];

      flux[find] += *fscale * flux0 * flux1;
    }
  }
}
