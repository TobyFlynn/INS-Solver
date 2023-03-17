inline void ins_3d_vis_x(const DG_FP *t, const DG_FP *g0, const int *bc_type, const int *faceNum,
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
      // bcs[i] = sin(PI * (*t));
      bcs[i] = 1.0;
    }
  } else if(*bc_type == LW_NO_SLIP_WALL_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      bcs[i] = 0.0;
    }
  } else if(*bc_type == LW_SLIP_WALL_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      const DG_FP dot = nx[0] * u[fmaskB[i]] + ny[0] * v[fmaskB[i]] + nz[0] * w[fmaskB[i]];
      bcs[i] = (u[fmaskB[i]] - dot * nx[0]) / *g0;
    }
  } else if(*bc_type == LW_OUTFLOW_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      bcs[i] = 0.0;
    }
  }
}
