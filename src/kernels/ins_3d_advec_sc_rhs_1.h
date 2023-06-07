inline void ins_3d_advec_sc_rhs_1(const int *faceNum, const int *fmaskL_corrected,
                           const int *fmaskR_corrected, const DG_FP *nx,
                           const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                           const DG_FP **us, const DG_FP **vs, const DG_FP **ws,
                           const DG_FP **ub, const DG_FP **vb, const DG_FP **wb,
                           DG_FP **f0, DG_FP **f1, DG_FP **f2) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];
  const int fIndL = faceNum[0] * DG_NPF;
  const int fIndR = faceNum[1] * DG_NPF;

  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];
    DG_FP lVel = nx[0] * ub[0][fmaskL_ind] + ny[0] * vb[0][fmaskL_ind] + nz[0] * wb[0][fmaskL_ind];
    DG_FP rVel = nx[1] * ub[1][fmaskR_ind] + ny[1] * vb[1][fmaskR_ind] + nz[1] * wb[1][fmaskR_ind];
    DG_FP vel = fmax(fabs(lVel), fabs(rVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Left numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];
    DG_FP f00L = ub[0][fmaskL_ind] * us[0][fmaskL_ind];
    DG_FP f01L = ub[0][fmaskL_ind] * vs[0][fmaskL_ind];
    DG_FP f02L = ub[0][fmaskL_ind] * ws[0][fmaskL_ind];
    DG_FP f10L = vb[0][fmaskL_ind] * us[0][fmaskL_ind];
    DG_FP f11L = vb[0][fmaskL_ind] * vs[0][fmaskL_ind];
    DG_FP f12L = vb[0][fmaskL_ind] * ws[0][fmaskL_ind];
    DG_FP f20L = wb[0][fmaskL_ind] * us[0][fmaskL_ind];
    DG_FP f21L = wb[0][fmaskL_ind] * vs[0][fmaskL_ind];
    DG_FP f22L = wb[0][fmaskL_ind] * ws[0][fmaskL_ind];

    DG_FP f00R = ub[1][fmaskR_ind] * us[1][fmaskR_ind];
    DG_FP f01R = ub[1][fmaskR_ind] * vs[1][fmaskR_ind];
    DG_FP f02R = ub[1][fmaskR_ind] * ws[1][fmaskR_ind];
    DG_FP f10R = vb[1][fmaskR_ind] * us[1][fmaskR_ind];
    DG_FP f11R = vb[1][fmaskR_ind] * vs[1][fmaskR_ind];
    DG_FP f12R = vb[1][fmaskR_ind] * ws[1][fmaskR_ind];
    DG_FP f20R = wb[1][fmaskR_ind] * us[1][fmaskR_ind];
    DG_FP f21R = wb[1][fmaskR_ind] * vs[1][fmaskR_ind];
    DG_FP f22R = wb[1][fmaskR_ind] * ws[1][fmaskR_ind];

    f0[0][fIndL + i] += 0.5 * fscale[0] * (-nx[0] * (f00L - f00R)
                        - ny[0] * (f01L - f01R) - nz[0] * (f02L - f02R)
                        - maxVel * (us[1][fmaskR_ind] - us[0][fmaskL_ind]));
    f1[0][fIndL + i] += 0.5 * fscale[0] * (-nx[0] * (f10L - f10R)
                        - ny[0] * (f11L - f11R) - nz[0] * (f12L - f12R)
                        - maxVel * (vs[1][fmaskR_ind] - vs[0][fmaskL_ind]));
    f2[0][fIndL + i] += 0.5 * fscale[0] * (-nx[0] * (f20L - f20R)
                        - ny[0] * (f21L - f21R) - nz[0] * (f22L - f22R)
                        - maxVel * (ws[1][fmaskR_ind] - ws[0][fmaskL_ind]));
  }

  // Right numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskR_ind = fmaskR[i];
    const int fmaskL_ind = fmaskL_corrected[i];
    DG_FP f00R = ub[1][fmaskR_ind] * us[1][fmaskR_ind];
    DG_FP f01R = ub[1][fmaskR_ind] * vs[1][fmaskR_ind];
    DG_FP f02R = ub[1][fmaskR_ind] * ws[1][fmaskR_ind];
    DG_FP f10R = vb[1][fmaskR_ind] * us[1][fmaskR_ind];
    DG_FP f11R = vb[1][fmaskR_ind] * vs[1][fmaskR_ind];
    DG_FP f12R = vb[1][fmaskR_ind] * ws[1][fmaskR_ind];
    DG_FP f20R = wb[1][fmaskR_ind] * us[1][fmaskR_ind];
    DG_FP f21R = wb[1][fmaskR_ind] * vs[1][fmaskR_ind];
    DG_FP f22R = wb[1][fmaskR_ind] * ws[1][fmaskR_ind];

    DG_FP f00L = ub[0][fmaskL_ind] * us[0][fmaskL_ind];
    DG_FP f01L = ub[0][fmaskL_ind] * vs[0][fmaskL_ind];
    DG_FP f02L = ub[0][fmaskL_ind] * ws[0][fmaskL_ind];
    DG_FP f10L = vb[0][fmaskL_ind] * us[0][fmaskL_ind];
    DG_FP f11L = vb[0][fmaskL_ind] * vs[0][fmaskL_ind];
    DG_FP f12L = vb[0][fmaskL_ind] * ws[0][fmaskL_ind];
    DG_FP f20L = wb[0][fmaskL_ind] * us[0][fmaskL_ind];
    DG_FP f21L = wb[0][fmaskL_ind] * vs[0][fmaskL_ind];
    DG_FP f22L = wb[0][fmaskL_ind] * ws[0][fmaskL_ind];

    f0[1][fIndR + i] += 0.5 * fscale[1] * (-nx[1] * (f00R - f00L)
                        - ny[1] * (f01R - f01L) - nz[1] * (f02R - f02L)
                        - maxVel * (us[0][fmaskL_ind] - us[1][fmaskR_ind]));
    f1[1][fIndR + i] += 0.5 * fscale[1] * (-nx[1] * (f10R - f10L)
                        - ny[1] * (f11R - f11L) - nz[1] * (f12R - f12L)
                        - maxVel * (vs[0][fmaskL_ind] - vs[1][fmaskR_ind]));
    f2[1][fIndR + i] += 0.5 * fscale[1] * (-nx[1] * (f20R - f20L)
                        - ny[1] * (f21R - f21L) - nz[1] * (f22R - f22L)
                        - maxVel * (ws[0][fmaskL_ind] - ws[1][fmaskR_ind]));
  }
}
