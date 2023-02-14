inline void ins_3d_advec_1(const int *faceNum, const int *fmaskL_corrected,
                           const int *fmaskR_corrected, const DG_FP *nx,
                           const DG_FP *ny, const DG_FP *nz,
                           const DG_FP *fscale, const DG_FP **u,
                           const DG_FP **v, const DG_FP **w, DG_FP **f0,
                           DG_FP **f1, DG_FP **f2) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];
  const int fIndL = faceNum[0] * DG_NPF;
  const int fIndR = faceNum[1] * DG_NPF;

  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    DG_FP lVel = nx[0] * u[0][fmaskL[i]] + ny[0] * v[0][fmaskL[i]] + nz[0] * w[0][fmaskL[i]];
    DG_FP rVel = nx[1] * u[1][fmaskR_corrected[i]] + ny[1] * v[1][fmaskR_corrected[i]] + nz[1] * w[1][fmaskR_corrected[i]];
    DG_FP vel = fmax(fabs(lVel), fabs(rVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Left numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    DG_FP f00L = u[0][fmaskL[i]] * u[0][fmaskL[i]];
    DG_FP f01L = u[0][fmaskL[i]] * v[0][fmaskL[i]];
    DG_FP f02L = u[0][fmaskL[i]] * w[0][fmaskL[i]];
    DG_FP f10L = v[0][fmaskL[i]] * u[0][fmaskL[i]];
    DG_FP f11L = v[0][fmaskL[i]] * v[0][fmaskL[i]];
    DG_FP f12L = v[0][fmaskL[i]] * w[0][fmaskL[i]];
    DG_FP f20L = w[0][fmaskL[i]] * u[0][fmaskL[i]];
    DG_FP f21L = w[0][fmaskL[i]] * v[0][fmaskL[i]];
    DG_FP f22L = w[0][fmaskL[i]] * w[0][fmaskL[i]];

    DG_FP f00R = u[1][fmaskR_corrected[i]] * u[1][fmaskR_corrected[i]];
    DG_FP f01R = u[1][fmaskR_corrected[i]] * v[1][fmaskR_corrected[i]];
    DG_FP f02R = u[1][fmaskR_corrected[i]] * w[1][fmaskR_corrected[i]];
    DG_FP f10R = v[1][fmaskR_corrected[i]] * u[1][fmaskR_corrected[i]];
    DG_FP f11R = v[1][fmaskR_corrected[i]] * v[1][fmaskR_corrected[i]];
    DG_FP f12R = v[1][fmaskR_corrected[i]] * w[1][fmaskR_corrected[i]];
    DG_FP f20R = w[1][fmaskR_corrected[i]] * u[1][fmaskR_corrected[i]];
    DG_FP f21R = w[1][fmaskR_corrected[i]] * v[1][fmaskR_corrected[i]];
    DG_FP f22R = w[1][fmaskR_corrected[i]] * w[1][fmaskR_corrected[i]];

    f0[0][fIndL + i] += 0.5 * fscale[0] * (-nx[0] * (f00L - f00R)
                        - ny[0] * (f01L - f01R) - nz[0] * (f02L - f02R)
                        - maxVel * (u[1][fmaskR_corrected[i]] - u[0][fmaskL[i]]));
    f1[0][fIndL + i] += 0.5 * fscale[0] * (-nx[0] * (f10L - f10R)
                        - ny[0] * (f11L - f11R) - nz[0] * (f12L - f12R)
                        - maxVel * (v[1][fmaskR_corrected[i]] - v[0][fmaskL[i]]));
    f2[0][fIndL + i] += 0.5 * fscale[0] * (-nx[0] * (f20L - f20R)
                        - ny[0] * (f21L - f21R) - nz[0] * (f22L - f22R)
                        - maxVel * (w[1][fmaskR_corrected[i]] - w[0][fmaskL[i]]));
  }

  // Right numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    DG_FP f00R = u[1][fmaskR[i]] * u[1][fmaskR[i]];
    DG_FP f01R = u[1][fmaskR[i]] * v[1][fmaskR[i]];
    DG_FP f02R = u[1][fmaskR[i]] * w[1][fmaskR[i]];
    DG_FP f10R = v[1][fmaskR[i]] * u[1][fmaskR[i]];
    DG_FP f11R = v[1][fmaskR[i]] * v[1][fmaskR[i]];
    DG_FP f12R = v[1][fmaskR[i]] * w[1][fmaskR[i]];
    DG_FP f20R = w[1][fmaskR[i]] * u[1][fmaskR[i]];
    DG_FP f21R = w[1][fmaskR[i]] * v[1][fmaskR[i]];
    DG_FP f22R = w[1][fmaskR[i]] * w[1][fmaskR[i]];

    DG_FP f00L = u[0][fmaskL_corrected[i]] * u[0][fmaskL_corrected[i]];
    DG_FP f01L = u[0][fmaskL_corrected[i]] * v[0][fmaskL_corrected[i]];
    DG_FP f02L = u[0][fmaskL_corrected[i]] * w[0][fmaskL_corrected[i]];
    DG_FP f10L = v[0][fmaskL_corrected[i]] * u[0][fmaskL_corrected[i]];
    DG_FP f11L = v[0][fmaskL_corrected[i]] * v[0][fmaskL_corrected[i]];
    DG_FP f12L = v[0][fmaskL_corrected[i]] * w[0][fmaskL_corrected[i]];
    DG_FP f20L = w[0][fmaskL_corrected[i]] * u[0][fmaskL_corrected[i]];
    DG_FP f21L = w[0][fmaskL_corrected[i]] * v[0][fmaskL_corrected[i]];
    DG_FP f22L = w[0][fmaskL_corrected[i]] * w[0][fmaskL_corrected[i]];

    f0[1][fIndR + i] += 0.5 * fscale[1] * (-nx[1] * (f00R - f00L)
                        - ny[1] * (f01R - f01L) - nz[1] * (f02R - f02L)
                        - maxVel * (u[0][fmaskL_corrected[i]] - u[1][fmaskR[i]]));
    f1[1][fIndR + i] += 0.5 * fscale[1] * (-nx[1] * (f10R - f10L)
                        - ny[1] * (f11R - f11L) - nz[1] * (f12R - f12L)
                        - maxVel * (v[0][fmaskL_corrected[i]] - v[1][fmaskR[i]]));
    f2[1][fIndR + i] += 0.5 * fscale[1] * (-nx[1] * (f20R - f20L)
                        - ny[1] * (f21R - f21L) - nz[1] * (f22R - f22L)
                        - maxVel * (w[0][fmaskL_corrected[i]] - w[1][fmaskR[i]]));
  }
}
