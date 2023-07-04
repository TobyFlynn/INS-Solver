inline void ins_3d_vis_y(const DG_FP *t, const DG_FP *g0, const int *bc_type, const int *faceNum,
                         const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                         const DG_FP *x, const DG_FP *y, const DG_FP *z,
                         const DG_FP *u, const DG_FP *v, const DG_FP *w,
                         int *type, DG_FP *bcs) {
  if(*bc_type == LW_OUTFLOW_BC)
    *type = 1;
  else
    *type = 0;

  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * DG_NPF];
  const DG_FP PI = 3.141592653589793238463;
  if(*bc_type == LW_INFLOW_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      bcs[i] = 0.0;
    }
  } else if(*bc_type == LW_NO_SLIP_WALL_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      bcs[i] = 0.0;
    }
  } else if(*bc_type == LW_SLIP_WALL_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      const DG_FP mag = sqrt(u[fmask_ind] * u[fmask_ind] + v[fmask_ind] * v[fmask_ind] + w[fmask_ind] * w[fmask_ind]);
      const DG_FP dot = nx[0] * u[fmask_ind] + ny[0] * v[fmask_ind] + nz[0] * w[fmask_ind];
      const DG_FP uR = u[fmask_ind] - dot * nx[0];
      const DG_FP vR = v[fmask_ind] - dot * ny[0];
      const DG_FP wR = w[fmask_ind] - dot * nz[0];
      const DG_FP mag2 = sqrt(uR * uR + vR * vR + wR * wR);
      const DG_FP mag_factor = fabs(mag) < 1e-8 || fabs(mag2) < 1e-8 ? 1.0 : mag / mag2;
      bcs[i] = (mag_factor * vR) / *g0;
    }
  } else if(*bc_type == LW_OUTFLOW_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      bcs[i] = 0.0;
    }
  }
}
