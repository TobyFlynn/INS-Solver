inline void ls_2d_stencil_avg_0(const int *faceNum, const bool *reverse,
                  const DG_FP **s, const DG_FP **stencil, DG_FP **avg) {
  const bool rev = *reverse;
  const int edgeL = faceNum[0];
  const int edgeR = faceNum[1];
  const int *fmaskL = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeL * DG_NPF];
  const int *fmaskR = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeR * DG_NPF];

  if(stencil[0][0] == stencil[1][0])
    return;

  if(stencil[0][0] == 1.0) {
    for(int i = 0; i < DG_NPF; i++) {
      const int indL = fmaskL[i];
      const int indR = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];
      avg[0][edgeL * DG_NPF + i] = 0.5 * (s[0][indL] + s[1][indR]);
    }
  }

  if(stencil[1][0] == 1.0) {
    for(int i = 0; i < DG_NPF; i++) {
      const int indR = fmaskR[i];
      const int indL = rev ? fmaskL[DG_NPF - i - 1] : fmaskL[i];
      avg[1][edgeR * DG_NPF + i] = 0.5 * (s[0][indL] + s[1][indR]);
    }
  }
}