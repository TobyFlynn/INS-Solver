inline void ins_2d_vis_bc_y(const DG_FP *t, const DG_FP *g0, const int *bedge_type,
                            const int *bedgeNum, const DG_FP *nx, const DG_FP *ny,
                            const DG_FP *x, const DG_FP *y, const DG_FP *u, const DG_FP *v,
                            DG_FP *out) {
  const int edge = *bedgeNum;
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];
  if(*bedge_type == BC_TYPE_NO_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      out[i] = 0.0;
    }
  } else if(*bedge_type == BC_TYPE_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      const DG_FP u_scaled = u[fmask_ind] / *g0;
      const DG_FP v_scaled = v[fmask_ind] / *g0;
      const DG_FP mag = sqrt(u_scaled * u_scaled + v_scaled * v_scaled);
      const DG_FP dot = nx[0] * u_scaled + ny[0] * v_scaled;
      const DG_FP u_new = u_scaled - dot * nx[0];
      const DG_FP v_new = v_scaled - dot * ny[0];
      const DG_FP mag2 = sqrt(u_new * u_new + v_new * v_new);
      const DG_FP mag_factor = fabs(mag) < 1e-8 || fabs(mag2) < 1e-8 ? 1.0 : mag / mag2;
      out[i] = v_new * mag_factor;
    }
  } else if(*bedge_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      out[i] = 0.0;
    }
  } else {
    int vis_type;
    ps2d_custom_bc_get_vis_type(*bedge_type, vis_type);
    if(vis_type == BC_DIRICHLET) {
      for(int i = 0; i < DG_NPF; i++) {
        const int fmask_ind = fmask[i];
        const DG_FP u_scaled = u[fmask_ind] / *g0;
        const DG_FP v_scaled = v[fmask_ind] / *g0;
        DG_FP u_bc = 0.0;
        DG_FP v_bc = 0.0;
        ps2d_custom_bc_get_vel(*bedge_type, *t, x[fmask_ind], y[fmask_ind], nx[0], ny[0], u_scaled, v_scaled, u_bc, v_bc);
        out[i] = v_bc;
      }
    } else {
      for(int i = 0; i < DG_NPF; i++) {
        const int fmask_ind = fmask[i];
        DG_FP u_bc_neumann = 0.0;
        DG_FP v_bc_neumann = 0.0;
        ps2d_custom_bc_get_vis_neumann(*bedge_type, *t, x[fmask_ind], y[fmask_ind], nx[0], ny[0], u_bc_neumann, v_bc_neumann);
        out[i] = v_bc_neumann;
      }
    }
  }
}
