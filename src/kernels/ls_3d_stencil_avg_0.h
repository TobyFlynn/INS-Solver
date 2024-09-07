inline void ls_3d_stencil_avg_0(const int *faceNum, const int *fmaskL_corrected,
                                const int *fmaskR_corrected, const DG_FP **s, 
                                const DG_FP **stencil, DG_FP **avg) {
  const int faceL = faceNum[0];
  const int faceR = faceNum[1];
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceL * DG_NPF];
  const int *fmaskR = &fmask[faceR * DG_NPF];

  if(stencil[0][0] == stencil[1][0])
    return;

  if(stencil[0][0] == 1.0) {
    for(int i = 0; i < DG_NPF; i++) {
      const int indL = fmaskL[i];
      const int indR = fmaskR_corrected[i];
      avg[0][faceL * DG_NPF + i] = 0.5 * (s[0][indL] + s[1][indR]);
    }
  }

  if(stencil[1][0] == 1.0) {
    for(int i = 0; i < DG_NPF; i++) {
      const int indR = fmaskR[i];
      const int indL = fmaskL_corrected[i];
      avg[1][faceR * DG_NPF + i] = 0.5 * (s[0][indL] + s[1][indR]);
    }
  }
}