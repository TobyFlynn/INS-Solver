inline void ins_2d_advec_oi_1(const int *faceNum, const bool *reverse,
                        const DG_FP **u, const DG_FP **v, DG_FP **uM,
                        DG_FP **vM, DG_FP **uP, DG_FP **vP) {
  const bool rev = *reverse;
  const int edgeL = faceNum[0];
  const int edgeR = faceNum[1];
  const int *fmaskL = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeL * DG_NPF];
  const int *fmaskR = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeR * DG_NPF];

  for(int i = 0; i < DG_NPF; i++) {
    const int indL = fmaskL[i];
    const int indR = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];
    
    uM[0][edgeL * DG_NPF + i] = u[0][indL];
    vM[0][edgeL * DG_NPF + i] = v[0][indL];
    uP[0][edgeL * DG_NPF + i] = u[1][indR];
    vP[0][edgeL * DG_NPF + i] = v[1][indR];
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int indR = fmaskR[i];
    const int indL = rev ? fmaskL[DG_NPF - i - 1] : fmaskL[i];
    
    uM[1][edgeR * DG_NPF + i] = u[1][indR];
    vM[1][edgeR * DG_NPF + i] = v[1][indR];
    uP[1][edgeR * DG_NPF + i] = u[0][indL];
    vP[1][edgeR * DG_NPF + i] = v[0][indL];
  }
}
