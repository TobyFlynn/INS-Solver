inline void ins_2d_proj_6(const int *faceNum, const bool *reverse, const DG_FP *sJ,
                          const DG_FP **pen, const DG_FP **u, const DG_FP **v,
                          DG_FP **f0, DG_FP **f1) {
  const bool rev = *reverse;
  const int edgeL = faceNum[0];
  const int edgeR = faceNum[1];
  const int *fmaskL = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeL * DG_NPF];
  const int *fmaskR = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeR * DG_NPF];
  const int fIndL = edgeL * DG_NPF;
  const int fIndR = edgeR * DG_NPF;

  const DG_FP pen_f = 0.5 * (pen[0][0] + pen[1][0]);

  // Left numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];

    f0[0][fIndL + i] = sJ[0] * pen_f * (u[0][fmaskL_ind] - u[1][fmaskR_ind]);
    f1[0][fIndL + i] = sJ[0] * pen_f * (v[0][fmaskL_ind] - v[1][fmaskR_ind]);
  }

  // Right numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskR_ind = fmaskR[i];
    const int fmaskL_ind = rev ? fmaskL[DG_NPF - i - 1] : fmaskL[i];

    f0[1][fIndR + i] = sJ[1] * pen_f * (u[1][fmaskR_ind] - u[0][fmaskL_ind]);
    f1[1][fIndR + i] = sJ[1] * pen_f * (v[1][fmaskR_ind] - v[0][fmaskL_ind]);
  }
}
