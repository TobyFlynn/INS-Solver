inline void ins_3d_advec_2(const DG_FP *t, const int *bc_type, const int *faceNum,
                           const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                           const DG_FP *fscale, const DG_FP *x, const DG_FP *y,
                           const DG_FP *z, const DG_FP *u, const DG_FP *v,
                           const DG_FP *w, DG_FP *f0, DG_FP *f1, DG_FP *f2) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;
  const DG_FP PI = 3.141592653589793238463;
  DG_FP uR[DG_NPF], vR[DG_NPF], wR[DG_NPF];
  if(*bc_type == LW_INFLOW_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      // uR[i] = sin(PI * (*t));
      DG_FP tmp = y[fmaskB[i]] * y[fmaskB[i]] + z[fmaskB[i]] * z[fmaskB[i]] - LW_INLET_RADIUS * LW_INLET_RADIUS;
      tmp = fabs(tmp);
      tmp = fmin(1.0, tmp / 0.1);
      uR[i] = tmp * 1.0;
      vR[i] = 0.0;
      wR[i] = 0.0;
    }
  } else if(*bc_type == LW_OUTFLOW_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      uR[i] = u[fmaskB[i]];
      vR[i] = v[fmaskB[i]];
      wR[i] = w[fmaskB[i]];
    }
  } else if(*bc_type == LW_SLIP_WALL_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      const DG_FP dot = nx[0] * u[fmaskB[i]] + ny[0] * v[fmaskB[i]] + nz[0] * w[fmaskB[i]];
      uR[i] = u[fmaskB[i]] - dot * nx[0];
      vR[i] = v[fmaskB[i]] - dot * ny[0];
      wR[i] = w[fmaskB[i]] - dot * nz[0];
    }
  } else if(*bc_type == LW_NO_SLIP_WALL_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      uR[i] = 0.0;
      vR[i] = 0.0;
      wR[i] = 0.0;
    }
  }

  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    DG_FP lVel = *nx * u[fmaskB[i]] + *ny * v[fmaskB[i]] + *nz * w[fmaskB[i]];
    DG_FP rVel = *nx * uR[i] + *ny * vR[i] + *nz * wR[i];
    DG_FP vel = fmax(fabs(lVel), fabs(rVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Left numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    DG_FP f00L = u[fmaskB[i]] * u[fmaskB[i]];
    DG_FP f01L = u[fmaskB[i]] * v[fmaskB[i]];
    DG_FP f02L = u[fmaskB[i]] * w[fmaskB[i]];
    DG_FP f10L = v[fmaskB[i]] * u[fmaskB[i]];
    DG_FP f11L = v[fmaskB[i]] * v[fmaskB[i]];
    DG_FP f12L = v[fmaskB[i]] * w[fmaskB[i]];
    DG_FP f20L = w[fmaskB[i]] * u[fmaskB[i]];
    DG_FP f21L = w[fmaskB[i]] * v[fmaskB[i]];
    DG_FP f22L = w[fmaskB[i]] * w[fmaskB[i]];

    DG_FP f00R = uR[i] * uR[i];
    DG_FP f01R = uR[i] * vR[i];
    DG_FP f02R = uR[i] * wR[i];
    DG_FP f10R = vR[i] * uR[i];
    DG_FP f11R = vR[i] * vR[i];
    DG_FP f12R = vR[i] * wR[i];
    DG_FP f20R = wR[i] * uR[i];
    DG_FP f21R = wR[i] * vR[i];
    DG_FP f22R = wR[i] * wR[i];

    f0[fInd + i] += 0.5 * *fscale * (-*nx * (f00L - f00R)
                        - *ny * (f01L - f01R) - *nz * (f02L - f02R)
                        - maxVel * (uR[i] - u[fmaskB[i]]));
    f1[fInd + i] += 0.5 * *fscale * (-*nx * (f10L - f10R)
                        - *ny * (f11L - f11R) - *nz * (f12L - f12R)
                        - maxVel * (vR[i] - v[fmaskB[i]]));
    f2[fInd + i] += 0.5 * *fscale * (-*nx * (f20L - f20R)
                        - *ny * (f21L - f21R) - *nz * (f22L - f22R)
                        - maxVel * (wR[i] - w[fmaskB[i]]));
  }
}
