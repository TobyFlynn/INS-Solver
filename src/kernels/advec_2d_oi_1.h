inline void advec_2d_oi_1(const int *faceNum, const bool *reverse, const DG_FP **u, 
                          const DG_FP **v, const DG_FP **val, DG_FP **uM, DG_FP **vM, 
                          DG_FP **valM, DG_FP **valP) {
  const bool rev = *reverse;
  const int edgeL = faceNum[0];
  const int edgeR = faceNum[1];
  const int *fmaskL = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeL * DG_NPF];
  const int *fmaskR = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeR * DG_NPF];

  for(int i = 0; i < DG_NPF; i++) {
    const int indL = fmaskL[i];
    const int indR = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];
    
    uM[0][edgeL * DG_NPF + i]   = u[0][indL];
    vM[0][edgeL * DG_NPF + i]   = v[0][indL];
    valM[0][edgeL * DG_NPF + i] = val[0][indL];
    valP[0][edgeL * DG_NPF + i] = val[1][indR];
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int indR = fmaskR[i];
    const int indL = rev ? fmaskL[DG_NPF - i - 1] : fmaskL[i];
    
    uM[1][edgeR * DG_NPF + i]   = u[1][indR];
    vM[1][edgeR * DG_NPF + i]   = v[1][indR];
    valM[1][edgeR * DG_NPF + i] = val[1][indR];
    valP[1][edgeR * DG_NPF + i] = val[0][indL];
  }
}
