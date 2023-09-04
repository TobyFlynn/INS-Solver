inline void ins_2d_st_2(const int *faceNum, const bool *reverse,
                        const DG_FP **s, DG_FP **sM, DG_FP **sP) {
  const bool rev = *reverse;
  const int edgeL = faceNum[0];
  const int edgeR = faceNum[1];
  const int *fmaskL = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeL * DG_NPF];
  const int *fmaskR = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeR * DG_NPF];

  for(int i = 0; i < DG_NPF; i++) {
    const int indL = fmaskL[i];
    const int indR = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];
    
    sM[0][edgeL * DG_NPF + i] = s[0][indL];
    sP[0][edgeL * DG_NPF + i] = s[1][indR];
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int indR = fmaskR[i];
    const int indL = rev ? fmaskL[DG_NPF - i - 1] : fmaskL[i];
    
    sM[1][edgeR * DG_NPF + i] = s[1][indR];
    sP[1][edgeR * DG_NPF + i] = s[0][indL];
  }
}
