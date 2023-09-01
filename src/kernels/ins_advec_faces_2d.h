inline void ins_advec_faces_2d(const int *faceNum, const bool *reverse,
                               const DG_FP *nx, const DG_FP *ny,
                               const DG_FP *fscale, const DG_FP **u,
                               const DG_FP **v, DG_FP **f0, DG_FP **f1) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];
  const bool rev = *reverse;

  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];
    DG_FP lVel = nx[0] * u[0][fmaskL_ind] + ny[0] * v[0][fmaskL_ind];
    DG_FP rVel = nx[1] * u[1][fmaskR_ind] + ny[1] * v[1][fmaskR_ind];
    DG_FP vel = fmax(fabs(lVel), fabs(rVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    DG_FP f00L = u[0][fmaskL_ind] * u[0][fmaskL_ind];
    DG_FP f01L = u[0][fmaskL_ind] * v[0][fmaskL_ind];
    DG_FP f10L = v[0][fmaskL_ind] * u[0][fmaskL_ind];
    DG_FP f11L = v[0][fmaskL_ind] * v[0][fmaskL_ind];

    const int fmaskR_ind = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];
    DG_FP f00R = u[1][fmaskR_ind] * u[1][fmaskR_ind];
    DG_FP f01R = u[1][fmaskR_ind] * v[1][fmaskR_ind];
    DG_FP f10R = v[1][fmaskR_ind] * u[1][fmaskR_ind];
    DG_FP f11R = v[1][fmaskR_ind] * v[1][fmaskR_ind];

    const int fIndL = faceNum[0] * DG_NPF + i;
    const int fIndR = rev ? faceNum[1] * DG_NPF + DG_NPF - i - 1 : faceNum[1] * DG_NPF + i;

    // f0[0][fIndL] = 0.5 * fscale[0] * (-nx[0] * (f00L - f00R)
    //                - ny[0] * (f01L - f01R) - maxVel * (u[1][fmaskR_ind] - u[0][fmaskL_ind]));
    // f1[0][fIndL] = 0.5 * fscale[0] * (-nx[0] * (f10L - f10R)
    //                - ny[0] * (f11L - f11R) - maxVel * (v[1][fmaskR_ind] - v[0][fmaskL_ind]));

    // f0[1][fIndR] = 0.5 * fscale[1] * (-nx[1] * (f00R - f00L)
    //                - ny[1] * (f01R - f01L) - maxVel * (u[0][fmaskL_ind] - u[1][fmaskR_ind]));
    // f1[1][fIndR] = 0.5 * fscale[1] * (-nx[1] * (f10R - f10L)
    //                - ny[1] * (f11R - f11L) - maxVel * (v[0][fmaskL_ind] - v[1][fmaskR_ind]));

    f0[0][fIndL] = 0.5 * fscale[0] * (nx[0] * (f00L + f00R)
                   + ny[0] * (f01L + f01R) + maxVel * (u[0][fmaskL_ind] - u[1][fmaskR_ind]));
    f1[0][fIndL] = 0.5 * fscale[0] * (nx[0] * (f10L + f10R)
                   + ny[0] * (f11L + f11R) + maxVel * (v[0][fmaskL_ind] - v[1][fmaskR_ind]));

    f0[1][fIndR] = 0.5 * fscale[1] * (nx[1] * (f00R + f00L)
                   + ny[1] * (f01R + f01L) + maxVel * (u[1][fmaskR_ind] - u[0][fmaskL_ind]));
    f1[1][fIndR] = 0.5 * fscale[1] * (nx[1] * (f10R + f10L)
                   + ny[1] * (f11R + f11L) + maxVel * (v[1][fmaskR_ind] - v[0][fmaskL_ind]));
  }
}
