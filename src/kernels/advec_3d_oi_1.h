inline void advec_3d_oi_1(const int *faceNum, const int *fmaskL_corrected,
                const int *fmaskR_corrected, const DG_FP **u, const DG_FP **v, 
                const DG_FP **w, const DG_FP **val, DG_FP **mU, DG_FP **mV, 
                DG_FP **mW, DG_FP **mVal, DG_FP **pVal) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];
  const int fIndL = faceNum[0] * DG_NPF;
  const int fIndR = faceNum[1] * DG_NPF;

  // Left
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];
    mU[0][fIndL + i]   = u[0][fmaskL_ind];
    mV[0][fIndL + i]   = v[0][fmaskL_ind];
    mW[0][fIndL + i]   = w[0][fmaskL_ind];
    mVal[0][fIndL + i] = val[0][fmaskL_ind];
    pVal[0][fIndL + i] = val[1][fmaskR_ind];
  }

  // Right
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskR_ind = fmaskR[i];
    const int fmaskL_ind = fmaskL_corrected[i];
    mU[1][fIndR + i]   = u[1][fmaskR_ind];
    mV[1][fIndR + i]   = v[1][fmaskR_ind];
    mW[1][fIndR + i]   = w[1][fmaskR_ind];
    mVal[1][fIndR + i] = val[1][fmaskR_ind];
    pVal[1][fIndR + i] = val[0][fmaskL_ind];
  }
}