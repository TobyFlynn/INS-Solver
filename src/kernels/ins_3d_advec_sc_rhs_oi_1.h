inline void ins_3d_advec_sc_rhs_oi_1(const int *faceNum, const int *fmaskL_corrected,
                           const int *fmaskR_corrected, const DG_FP **us, const DG_FP **vs, 
                           const DG_FP **ws, const DG_FP **ub, const DG_FP **vb, const DG_FP **wb,
                           DG_FP **mUs, DG_FP **mVs, DG_FP **mWs, DG_FP **pUs, DG_FP **pVs, DG_FP **pWs,
                           DG_FP **mUb, DG_FP **mVb, DG_FP **mWb, DG_FP **pUb, DG_FP **pVb, DG_FP **pWb) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];
  const int fIndL = faceNum[0] * DG_NPF;
  const int fIndR = faceNum[1] * DG_NPF;

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];
    mUs[0][fIndL + i] = us[0][fmaskL_ind];
    mVs[0][fIndL + i] = vs[0][fmaskL_ind];
    mWs[0][fIndL + i] = ws[0][fmaskL_ind];
    pUs[0][fIndL + i] = us[1][fmaskR_ind];
    pVs[0][fIndL + i] = vs[1][fmaskR_ind];
    pWs[0][fIndL + i] = ws[1][fmaskR_ind];
    mUb[0][fIndL + i] = ub[0][fmaskL_ind];
    mVb[0][fIndL + i] = vb[0][fmaskL_ind];
    mWb[0][fIndL + i] = wb[0][fmaskL_ind];
    pUb[0][fIndL + i] = ub[1][fmaskR_ind];
    pVb[0][fIndL + i] = vb[1][fmaskR_ind];
    pWb[0][fIndL + i] = wb[1][fmaskR_ind];
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskR_ind = fmaskR[i];
    const int fmaskL_ind = fmaskL_corrected[i];
    mUs[1][fIndR + i] = us[1][fmaskR_ind];
    mVs[1][fIndR + i] = vs[1][fmaskR_ind];
    mWs[1][fIndR + i] = ws[1][fmaskR_ind];
    pUs[1][fIndR + i] = us[0][fmaskL_ind];
    pVs[1][fIndR + i] = vs[0][fmaskL_ind];
    pWs[1][fIndR + i] = ws[0][fmaskL_ind];
    mUb[1][fIndR + i] = ub[1][fmaskR_ind];
    mVb[1][fIndR + i] = vb[1][fmaskR_ind];
    mWb[1][fIndR + i] = wb[1][fmaskR_ind];
    pUb[1][fIndR + i] = ub[0][fmaskL_ind];
    pVb[1][fIndR + i] = vb[0][fmaskL_ind];
    pWb[1][fIndR + i] = wb[0][fmaskL_ind];
  }
}
