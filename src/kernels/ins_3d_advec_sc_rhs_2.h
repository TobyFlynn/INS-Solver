inline void ins_3d_advec_sc_rhs_2(const int *bc_type, const int *faceNum,
                           const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                           const DG_FP *fscale, const DG_FP *x, const DG_FP *y,
                           const DG_FP *z, const DG_FP *us, const DG_FP *vs,
                           const DG_FP *ws, const DG_FP *ub, const DG_FP *vb,
                           const DG_FP *wb, DG_FP *f0, DG_FP *f1, DG_FP *f2) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;
  const DG_FP PI = 3.141592653589793238463;
  DG_FP usR[DG_NPF], vsR[DG_NPF], wsR[DG_NPF];
  DG_FP ubR[DG_NPF], vbR[DG_NPF], wbR[DG_NPF];
  if(*bc_type == LW_INFLOW_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      // uR[i] = sin(PI * (*t));
      DG_FP tmp = y[fmaskB[i]] * y[fmaskB[i]] + z[fmaskB[i]] * z[fmaskB[i]] - LW_INLET_RADIUS * LW_INLET_RADIUS;
      tmp = fabs(tmp);
      tmp = fmin(1.0, tmp / 0.1);
      usR[i] = tmp * 1.0;
      vsR[i] = 0.0;
      wsR[i] = 0.0;
      ubR[i] = tmp * 1.0;
      vbR[i] = 0.0;
      wbR[i] = 0.0;
    }
  } else if(*bc_type == LW_OUTFLOW_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      usR[i] = us[fmaskB[i]];
      vsR[i] = vs[fmaskB[i]];
      wsR[i] = ws[fmaskB[i]];
      ubR[i] = ub[fmaskB[i]];
      vbR[i] = vb[fmaskB[i]];
      wbR[i] = wb[fmaskB[i]];
    }
  } else if(*bc_type == LW_SLIP_WALL_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      const DG_FP dots = nx[0] * us[fmaskB[i]] + ny[0] * vs[fmaskB[i]] + nz[0] * ws[fmaskB[i]];
      usR[i] = us[fmaskB[i]] - dots * nx[0];
      vsR[i] = vs[fmaskB[i]] - dots * ny[0];
      wsR[i] = ws[fmaskB[i]] - dots * nz[0];
      const DG_FP dotb = nx[0] * ub[fmaskB[i]] + ny[0] * vb[fmaskB[i]] + nz[0] * wb[fmaskB[i]];
      ubR[i] = ub[fmaskB[i]] - dotb * nx[0];
      vbR[i] = vb[fmaskB[i]] - dotb * ny[0];
      wbR[i] = wb[fmaskB[i]] - dotb * nz[0];
    }
  } else if(*bc_type == LW_NO_SLIP_WALL_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      usR[i] = 0.0;
      vsR[i] = 0.0;
      wsR[i] = 0.0;
      ubR[i] = 0.0;
      vbR[i] = 0.0;
      wbR[i] = 0.0;
    }
  }

  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    DG_FP lVel = *nx * ub[fmaskB[i]] + *ny * vb[fmaskB[i]] + *nz * wb[fmaskB[i]];
    DG_FP rVel = *nx * ubR[i] + *ny * vbR[i] + *nz * wbR[i];
    DG_FP vel = fmax(fabs(lVel), fabs(rVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Left numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    DG_FP f00L = ub[fmaskB[i]] * us[fmaskB[i]];
    DG_FP f01L = ub[fmaskB[i]] * vs[fmaskB[i]];
    DG_FP f02L = ub[fmaskB[i]] * ws[fmaskB[i]];
    DG_FP f10L = vb[fmaskB[i]] * us[fmaskB[i]];
    DG_FP f11L = vb[fmaskB[i]] * vs[fmaskB[i]];
    DG_FP f12L = vb[fmaskB[i]] * ws[fmaskB[i]];
    DG_FP f20L = wb[fmaskB[i]] * us[fmaskB[i]];
    DG_FP f21L = wb[fmaskB[i]] * vs[fmaskB[i]];
    DG_FP f22L = wb[fmaskB[i]] * ws[fmaskB[i]];

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
                        - maxVel * (usR[i] - us[fmaskB[i]]));
    f1[fInd + i] += 0.5 * *fscale * (-*nx * (f10L - f10R)
                        - *ny * (f11L - f11R) - *nz * (f12L - f12R)
                        - maxVel * (vsR[i] - vs[fmaskB[i]]));
    f2[fInd + i] += 0.5 * *fscale * (-*nx * (f20L - f20R)
                        - *ny * (f21L - f21R) - *nz * (f22L - f22R)
                        - maxVel * (wsR[i] - ws[fmaskB[i]]));
  }
}
