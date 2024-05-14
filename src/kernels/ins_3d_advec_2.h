inline void ins_3d_advec_2(const DG_FP *t, const int *bc_type, const int *faceNum,
                           const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                           const DG_FP *fscale, const DG_FP *x, const DG_FP *y,
                           const DG_FP *z, const DG_FP *u, const DG_FP *v,
                           const DG_FP *w, DG_FP *f0, DG_FP *f1, DG_FP *f2) {
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;
  DG_FP uR[DG_NPF], vR[DG_NPF], wR[DG_NPF];
  if(*bc_type == BC_TYPE_NO_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      uR[i] = 0.0;
      vR[i] = 0.0;
      wR[i] = 0.0;
    }
  } else if(*bc_type == BC_TYPE_SLIP) {
    const DG_FP mag_normal = sqrt(*nx * *nx + *ny * *ny + *nz * *nz);
    const DG_FP nx_ = *nx / mag_normal;
    const DG_FP ny_ = *ny / mag_normal;
    const DG_FP nz_ = *nz / mag_normal;
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      const DG_FP dot = nx_ * u[fmask_ind] + ny_ * v[fmask_ind] + nz_ * w[fmask_ind];
      uR[i] = u[fmask_ind] - dot * nx_;
      vR[i] = v[fmask_ind] - dot * ny_;
      wR[i] = w[fmask_ind] - dot * nz_;
    }
  } else if(*bc_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      uR[i] = u[fmask_ind];
      vR[i] = v[fmask_ind];
      wR[i] = w[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      ps3d_custom_bc_get_vel(*bc_type, *t, x[fmask_ind], y[fmask_ind], z[fmask_ind],
                             *nx, *ny, *nz, u[fmask_ind], v[fmask_ind], w[fmask_ind],
                             uR[i], vR[i], wR[i]);
    }
  }

  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    DG_FP lVel = *nx * u[fmask_ind] + *ny * v[fmask_ind] + *nz * w[fmask_ind];
    DG_FP rVel = *nx * uR[i] + *ny * vR[i] + *nz * wR[i];
    DG_FP vel = fmax(fabs(lVel), fabs(rVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Left numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    DG_FP f00L = u[fmask_ind] * u[fmask_ind];
    DG_FP f01L = u[fmask_ind] * v[fmask_ind];
    DG_FP f02L = u[fmask_ind] * w[fmask_ind];
    DG_FP f10L = v[fmask_ind] * u[fmask_ind];
    DG_FP f11L = v[fmask_ind] * v[fmask_ind];
    DG_FP f12L = v[fmask_ind] * w[fmask_ind];
    DG_FP f20L = w[fmask_ind] * u[fmask_ind];
    DG_FP f21L = w[fmask_ind] * v[fmask_ind];
    DG_FP f22L = w[fmask_ind] * w[fmask_ind];

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
                        - maxVel * (uR[i] - u[fmask_ind]));
    f1[fInd + i] += 0.5 * *fscale * (-*nx * (f10L - f10R)
                        - *ny * (f11L - f11R) - *nz * (f12L - f12R)
                        - maxVel * (vR[i] - v[fmask_ind]));
    f2[fInd + i] += 0.5 * *fscale * (-*nx * (f20L - f20R)
                        - *ny * (f21L - f21R) - *nz * (f22L - f22R)
                        - maxVel * (wR[i] - w[fmask_ind]));
  }
}
