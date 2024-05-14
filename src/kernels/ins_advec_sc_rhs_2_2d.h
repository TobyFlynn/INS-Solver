inline void ins_advec_sc_rhs_2_2d(const DG_FP *t, const int *bedge_type, const int *bedgeNum,
                            const DG_FP *nx, const DG_FP *ny, const DG_FP *fscale,
                            const DG_FP *x, const DG_FP *y, const DG_FP *us,
                            const DG_FP *vs, const DG_FP *ub, const DG_FP *vb,
                            DG_FP *flux0, DG_FP *flux1) {
  DG_FP pUs[DG_NPF], pVs[DG_NPF], pUb[DG_NPF], pVb[DG_NPF];
  const int edgeNum = bedgeNum[0];
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeNum * DG_NPF];

  // Set boundary velocities
  if(*bedge_type == BC_TYPE_NO_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      pUs[i] = 0.0;
      pVs[i] = 0.0;
      pUb[i] = 0.0;
      pVb[i] = 0.0;
    }
  } else if (*bedge_type == BC_TYPE_SLIP) {
    DG_FP tangent_x = *ny;
    DG_FP tangent_y = -*nx;
    DG_FP tangent_mag = sqrt(tangent_x * tangent_x + tangent_y * tangent_y);
    tangent_x /= tangent_mag;
    tangent_y /= tangent_mag;
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      const DG_FP dot_s = us[fmask_ind] * tangent_x + vs[fmask_ind] * tangent_y;
      const DG_FP dot_b = ub[fmask_ind] * tangent_x + vb[fmask_ind] * tangent_y;
      pUs[i] = dot_s * tangent_x;
      pVs[i] = dot_s * tangent_y;
      pUb[i] = dot_b * tangent_x;
      pVb[i] = dot_b * tangent_y;
    }
  } else if(*bedge_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      pUs[i] = us[fmask_ind];
      pVs[i] = vs[fmask_ind];
      pUb[i] = ub[fmask_ind];
      pVb[i] = vb[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      ps2d_custom_bc_get_vel(*bedge_type, *t, x[fmask_ind], y[fmask_ind], nx[0], ny[0], us[fmask_ind], vs[fmask_ind], pUs[i], pVs[i]);
      ps2d_custom_bc_get_vel(*bedge_type, *t, x[fmask_ind], y[fmask_ind], nx[0], ny[0], ub[fmask_ind], vb[fmask_ind], pUb[i], pVb[i]);
    }
  }

  // Calculate numerical flux
  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];

    DG_FP mVel = ub[fmask_ind] * nx[0] + vb[fmask_ind] * ny[0];
    DG_FP pVel = pUb[i] * -nx[0] + pVb[i] * -ny[0];
    DG_FP vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    // Get interior flux terms
    DG_FP mF0 = ub[fmask_ind] * us[fmask_ind];
    DG_FP mF1 = ub[fmask_ind] * vs[fmask_ind];
    DG_FP mF2 = ub[fmask_ind] * vs[fmask_ind];
    DG_FP mF3 = vb[fmask_ind] * vs[fmask_ind];
    // Get exterior flux terms
    DG_FP pF0 = pUb[i] * pUs[i];
    DG_FP pF1 = pUb[i] * pVs[i];
    DG_FP pF2 = pUb[i] * pVs[i];
    DG_FP pF3 = pVb[i] * pVs[i];
    // Numerical flux
    const int flux_ind = edgeNum * DG_NPF + i;
    flux0[flux_ind] = 0.5 * fscale[0] * (-nx[0] * (mF0 - pF0) - ny[0] * (mF1 - pF1) - maxVel * (pUs[i] - us[fmask_ind]));
    flux1[flux_ind] = 0.5 * fscale[0] * (-nx[0] * (mF2 - pF2) - ny[0] * (mF3 - pF3) - maxVel * (pVs[i] - vs[fmask_ind]));
  }
}
