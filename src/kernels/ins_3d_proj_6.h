inline void ins_3d_proj_6(const int *faceNum, const int *fmaskL_corrected,
                          const int *fmaskR_corrected, const DG_FP *sJ,
                          const DG_FP **pen, const DG_FP **u, const DG_FP **v,
                          const DG_FP **w, DG_FP **f0, DG_FP **f1, DG_FP **f2) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];
  const int fIndL = faceNum[0] * DG_NPF;
  const int fIndR = faceNum[1] * DG_NPF;

  const DG_FP pen_f = 0.5 * (pen[0][0] + pen[1][0]);

  // Left numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];

    f0[0][fIndL + i] = sJ[0] * pen_f * (u[0][fmaskL_ind] - u[1][fmaskR_ind]);
    f1[0][fIndL + i] = sJ[0] * pen_f * (v[0][fmaskL_ind] - v[1][fmaskR_ind]);
    f2[0][fIndL + i] = sJ[0] * pen_f * (w[0][fmaskL_ind] - w[1][fmaskR_ind]);
  }

  // Right numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskR_ind = fmaskR[i];
    const int fmaskL_ind = fmaskL_corrected[i];

    f0[1][fIndR + i] = sJ[1] * pen_f * (u[1][fmaskR_ind] - u[0][fmaskL_ind]);
    f1[1][fIndR + i] = sJ[1] * pen_f * (v[1][fmaskR_ind] - v[0][fmaskL_ind]);
    f2[1][fIndR + i] = sJ[1] * pen_f * (w[1][fmaskR_ind] - w[0][fmaskL_ind]);
  }
}
