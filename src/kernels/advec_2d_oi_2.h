inline void advec_2d_oi_2(const int *bedgeNum, const DG_FP *u, const DG_FP *v,
                          const DG_FP *val, DG_FP *uM, DG_FP *vM, DG_FP *valM, 
                          DG_FP *valP) {
  const int edge = bedgeNum[0];
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];

  // Simple for surface tension case
  for(int i = 0; i < DG_NPF; i++) {
    const int find = fmask[i];
    
    uM[edge * DG_NPF + i]   = u[find];
    vM[edge * DG_NPF + i]   = v[find];
    valM[edge * DG_NPF + i] = val[find];
    valP[edge * DG_NPF + i] = val[find];
  }
}
