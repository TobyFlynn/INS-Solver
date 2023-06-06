inline void ins_advec_sc_rhs_1_2d(const int *faceNum, const bool *reverse,
                               const DG_FP *nx, const DG_FP *ny,
                               const DG_FP *fscale, const DG_FP **us,
                               const DG_FP **vs, const DG_FP **ub,
                               const DG_FP **vb, DG_FP **f0, DG_FP **f1) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];
  const bool rev = *reverse;

  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];
    DG_FP lVel = nx[0] * ub[0][fmaskL_ind] + ny[0] * vb[0][fmaskL_ind];
    DG_FP rVel = nx[1] * ub[1][fmaskR_ind] + ny[1] * vb[1][fmaskR_ind];
    DG_FP vel = fmax(fabs(lVel), fabs(rVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    DG_FP f00L = ub[0][fmaskL_ind] * us[0][fmaskL_ind];
    DG_FP f01L = ub[0][fmaskL_ind] * vs[0][fmaskL_ind];
    DG_FP f10L = vb[0][fmaskL_ind] * us[0][fmaskL_ind];
    DG_FP f11L = vb[0][fmaskL_ind] * vs[0][fmaskL_ind];

    const int fmaskR_ind = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];
    DG_FP f00R = ub[1][fmaskR_ind] * us[1][fmaskR_ind];
    DG_FP f01R = ub[1][fmaskR_ind] * vs[1][fmaskR_ind];
    DG_FP f10R = vb[1][fmaskR_ind] * us[1][fmaskR_ind];
    DG_FP f11R = vb[1][fmaskR_ind] * vs[1][fmaskR_ind];

    const int fIndL = faceNum[0] * DG_NPF + i;
    const int fIndR = rev ? faceNum[1] * DG_NPF + DG_NPF - i - 1 : faceNum[1] * DG_NPF + i;

    f0[0][fIndL] = 0.5 * fscale[0] * (-nx[0] * (f00L - f00R)
                   - ny[0] * (f01L - f01R) - maxVel * (us[1][fmaskR_ind] - us[0][fmaskL_ind]));
    f1[0][fIndL] = 0.5 * fscale[0] * (-nx[0] * (f10L - f10R)
                   - ny[0] * (f11L - f11R) - maxVel * (vs[1][fmaskR_ind] - vs[0][fmaskL_ind]));

    f0[1][fIndR] = 0.5 * fscale[1] * (-nx[1] * (f00R - f00L)
                   - ny[1] * (f01R - f01L) - maxVel * (us[0][fmaskL_ind] - us[1][fmaskR_ind]));
    f1[1][fIndR] = 0.5 * fscale[1] * (-nx[1] * (f10R - f10L)
                   - ny[1] * (f11R - f11L) - maxVel * (vs[0][fmaskL_ind] - vs[1][fmaskR_ind]));
  }
}
