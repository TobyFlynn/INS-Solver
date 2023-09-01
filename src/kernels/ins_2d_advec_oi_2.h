inline void ins_2d_advec_oi_2(const int *bedge_type, const int *bedgeNum,
                              const DG_FP *u, const DG_FP *v, DG_FP *uM, 
                              DG_FP *vM, DG_FP *uP, DG_FP *vP) {
  const int edge = bedgeNum[0];
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];

  // Simple for surface tension case
  for(int i = 0; i < DG_NPF; i++) {
    const int find = fmask[i];
    
    uM[edge * DG_NPF + i] = u[find];
    vM[edge * DG_NPF + i] = v[find];
    uP[edge * DG_NPF + i] = u[find];
    vP[edge * DG_NPF + i] = v[find];
  }
}
