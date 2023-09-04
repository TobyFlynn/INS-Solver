inline void ins_2d_st_3(const int *bedge_type, const int *bedgeNum,
                        const DG_FP *s, DG_FP *sM, DG_FP *sP) {
  const int edge = bedgeNum[0];
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];

  // Simple for surface tension case
  for(int i = 0; i < DG_NPF; i++) {
    const int find = fmask[i];
    
    sM[edge * DG_NPF + i] = s[find];
    sP[edge * DG_NPF + i] = s[find];
  }
}
