inline void ins_3d_advec_2(const double *t, const int *bc_type, const int *faceNum,
                           const double *nx, const double *ny, const double *nz,
                           const double *fscale, const double *x, const double *y,
                           const double *z, const double *u, const double *v,
                           const double *w, double *f0, double *f1, double *f2) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;
  const double PI = 3.141592653589793238463;
  double uR[DG_NPF], vR[DG_NPF], wR[DG_NPF];
  if(*bc_type == 0) {
    for(int i = 0; i < DG_NPF; i++) {
      // uR[i] = sin(PI * (*t)) * (y[fmask[i]] * (1.0 - y[fmask[i]]))  * (z[fmask[i]] * (1.0 - z[fmask[i]]));
      uR[i] = sin(PI * (*t));
      vR[i] = 0.0;
      wR[i] = 0.0;
    }
  } else if(*bc_type == 1) {
    for(int i = 0; i < DG_NPF; i++) {
      uR[i] = u[fmaskB[i]];
      vR[i] = v[fmaskB[i]];
      wR[i] = w[fmaskB[i]];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      uR[i] = 0.0;
      vR[i] = 0.0;
      wR[i] = 0.0;
    }
  }

  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    double lVel = *nx * u[fmaskB[i]] + *ny * v[fmaskB[i]] + *nz * w[fmaskB[i]];
    double rVel = *nx * uR[i] + *ny * vR[i] + *nz * wR[i];
    double vel = fmax(fabs(lVel), fabs(rVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Left numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    double f00L = u[fmaskB[i]] * u[fmaskB[i]];
    double f01L = u[fmaskB[i]] * v[fmaskB[i]];
    double f02L = u[fmaskB[i]] * w[fmaskB[i]];
    double f10L = v[fmaskB[i]] * u[fmaskB[i]];
    double f11L = v[fmaskB[i]] * v[fmaskB[i]];
    double f12L = v[fmaskB[i]] * w[fmaskB[i]];
    double f20L = w[fmaskB[i]] * u[fmaskB[i]];
    double f21L = w[fmaskB[i]] * v[fmaskB[i]];
    double f22L = w[fmaskB[i]] * w[fmaskB[i]];

    double f00R = uR[i] * uR[i];
    double f01R = uR[i] * vR[i];
    double f02R = uR[i] * wR[i];
    double f10R = vR[i] * uR[i];
    double f11R = vR[i] * vR[i];
    double f12R = vR[i] * wR[i];
    double f20R = wR[i] * uR[i];
    double f21R = wR[i] * vR[i];
    double f22R = wR[i] * wR[i];

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
