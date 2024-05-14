inline void ins_3d_advec_sc_rhs_2(const DG_FP *t, const int *bc_type, const int *faceNum,
                           const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                           const DG_FP *fscale, const DG_FP *x, const DG_FP *y,
                           const DG_FP *z, const DG_FP *us, const DG_FP *vs,
                           const DG_FP *ws, const DG_FP *ub, const DG_FP *vb,
                           const DG_FP *wb, DG_FP *f0, DG_FP *f1, DG_FP *f2) {
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;

  DG_FP usR[DG_NPF], vsR[DG_NPF], wsR[DG_NPF];
  DG_FP ubR[DG_NPF], vbR[DG_NPF], wbR[DG_NPF];
  if(*bc_type == BC_TYPE_NO_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      usR[i] = 0.0;
      vsR[i] = 0.0;
      wsR[i] = 0.0;
      ubR[i] = 0.0;
      vbR[i] = 0.0;
      wbR[i] = 0.0;
    }
  } else if(*bc_type == BC_TYPE_SLIP) {
    const DG_FP mag_normal = sqrt(*nx * *nx + *ny * *ny + *nz * *nz);
    const DG_FP nx_ = *nx / mag_normal;
    const DG_FP ny_ = *ny / mag_normal;
    const DG_FP nz_ = *nz / mag_normal;
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      // S
      const DG_FP dot_s = nx_ * us[fmask_ind] + ny_ * vs[fmask_ind] + nz_ * ws[fmask_ind];
      usR[i] = us[fmask_ind] - dot_s * nx_;
      vsR[i] = vs[fmask_ind] - dot_s * ny_;
      wsR[i] = ws[fmask_ind] - dot_s * nz_;
      // B
      const DG_FP dot_b = nx_ * ub[fmask_ind] + ny_ * vb[fmask_ind] + nz_ * wb[fmask_ind];
      ubR[i] = ub[fmask_ind] - dot_b * nx_;
      vbR[i] = vb[fmask_ind] - dot_b * ny_;
      wbR[i] = wb[fmask_ind] - dot_b * nz_;
    }
  } else if(*bc_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      usR[i] = us[fmask_ind];
      vsR[i] = vs[fmask_ind];
      wsR[i] = ws[fmask_ind];
      ubR[i] = ub[fmask_ind];
      vbR[i] = vb[fmask_ind];
      wbR[i] = wb[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      ps3d_custom_bc_get_vel(*bc_type, *t, x[fmask_ind], y[fmask_ind], z[fmask_ind],
                             *nx, *ny, *nz, us[fmask_ind], vs[fmask_ind], ws[fmask_ind],
                             usR[i], vsR[i], wsR[i]);
      ps3d_custom_bc_get_vel(*bc_type, *t, x[fmask_ind], y[fmask_ind], z[fmask_ind],
                             *nx, *ny, *nz, ub[fmask_ind], vb[fmask_ind], wb[fmask_ind],
                             ubR[i], vbR[i], wbR[i]);
    }
  }

  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    DG_FP lVel = *nx * ub[fmask_ind] + *ny * vb[fmask_ind] + *nz * wb[fmask_ind];
    DG_FP rVel = *nx * ubR[i] + *ny * vbR[i] + *nz * wbR[i];
    DG_FP vel = fmax(fabs(lVel), fabs(rVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Left numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    DG_FP f00L = ub[fmask_ind] * us[fmask_ind];
    DG_FP f01L = ub[fmask_ind] * vs[fmask_ind];
    DG_FP f02L = ub[fmask_ind] * ws[fmask_ind];
    DG_FP f10L = vb[fmask_ind] * us[fmask_ind];
    DG_FP f11L = vb[fmask_ind] * vs[fmask_ind];
    DG_FP f12L = vb[fmask_ind] * ws[fmask_ind];
    DG_FP f20L = wb[fmask_ind] * us[fmask_ind];
    DG_FP f21L = wb[fmask_ind] * vs[fmask_ind];
    DG_FP f22L = wb[fmask_ind] * ws[fmask_ind];

    DG_FP f00R = ubR[i] * usR[i];
    DG_FP f01R = ubR[i] * vsR[i];
    DG_FP f02R = ubR[i] * wsR[i];
    DG_FP f10R = vbR[i] * usR[i];
    DG_FP f11R = vbR[i] * vsR[i];
    DG_FP f12R = vbR[i] * wsR[i];
    DG_FP f20R = wbR[i] * usR[i];
    DG_FP f21R = wbR[i] * vsR[i];
    DG_FP f22R = wbR[i] * wsR[i];

    f0[fInd + i] += 0.5 * *fscale * (-*nx * (f00L - f00R)
                        - *ny * (f01L - f01R) - *nz * (f02L - f02R)
                        - maxVel * (usR[i] - us[fmask_ind]));
    f1[fInd + i] += 0.5 * *fscale * (-*nx * (f10L - f10R)
                        - *ny * (f11L - f11R) - *nz * (f12L - f12R)
                        - maxVel * (vsR[i] - vs[fmask_ind]));
    f2[fInd + i] += 0.5 * *fscale * (-*nx * (f20L - f20R)
                        - *ny * (f21L - f21R) - *nz * (f22L - f22R)
                        - maxVel * (wsR[i] - ws[fmask_ind]));
  }
}
