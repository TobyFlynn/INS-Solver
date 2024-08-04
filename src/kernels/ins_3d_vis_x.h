inline void ins_3d_vis_x(const DG_FP *t, const DG_FP *g0, const int *bc_type, const int *faceNum,
                         const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                         const DG_FP *x, const DG_FP *y, const DG_FP *z,
                         const DG_FP *u, const DG_FP *v, const DG_FP *w,
                         int *type, DG_FP *bcs) {
  if(*bc_type == BC_TYPE_NO_SLIP) {
    *type = BC_DIRICHLET;
  } else if(*bc_type == BC_TYPE_SLIP) {
    *type = BC_SLIP;
  } else if(*bc_type == BC_TYPE_NATURAL_OUTFLOW) {
    *type = BC_NEUMANN;
  } else {
    ps3d_custom_bc_get_vis_type(*bc_type, *type);
  }

  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  if(*bc_type == BC_TYPE_NO_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      bcs[i] = 0.0;
    }
  } else if(*bc_type == BC_TYPE_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      bcs[i] = 0.0;
    }
  } else if(*bc_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      bcs[i] = 0.0;
    }
  } else {
    if(*type == BC_DIRICHLET) {
      for(int i = 0; i < DG_NPF; i++) {
        const int fmask_ind = fmask[i];
        const DG_FP u_scaled = u[fmask_ind] / *g0;
        const DG_FP v_scaled = v[fmask_ind] / *g0;
        const DG_FP w_scaled = w[fmask_ind] / *g0;
        DG_FP u_new = 0.0;
        DG_FP v_new = 0.0;
        DG_FP w_new = 0.0;
        ps3d_custom_bc_get_vel(*bc_type, *t, x[fmask_ind], y[fmask_ind], z[fmask_ind],
                               *nx, *ny, *nz, u_scaled, v_scaled, w_scaled,
                               u_new, v_new, w_new);
        bcs[i] = u_new;
      }
    } else {
      for(int i = 0; i < DG_NPF; i++) {
        const int fmask_ind = fmask[i];
        const DG_FP u_scaled = u[fmask_ind] / *g0;
        const DG_FP v_scaled = v[fmask_ind] / *g0;
        const DG_FP w_scaled = w[fmask_ind] / *g0;
        DG_FP u_new = 0.0;
        DG_FP v_new = 0.0;
        DG_FP w_new = 0.0;
        ps3d_custom_bc_get_vis_neumann(*bc_type, *t, x[fmask_ind], y[fmask_ind], z[fmask_ind],
                               *nx, *ny, *nz, u_scaled, v_scaled, w_scaled,
                               u_new, v_new, w_new);
        bcs[i] = u_new;
      }
    }
  }
}
