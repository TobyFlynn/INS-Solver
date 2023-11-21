inline void ins_2d_st_3(const int *bedge_type, const int *bedgeNum,
                        const DG_FP *x, const DG_FP *y, const DG_FP *s,
                        DG_FP *sM, DG_FP *sP) {
  const int edge = bedgeNum[0];
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];

  // Simple for surface tension case
  for(int i = 0; i < DG_NPF; i++) {
    const int find = fmask[i];
    sM[edge * DG_NPF + i] = s[find];
  }

  if(*bedge_type == BC_TYPE_NO_SLIP || *bedge_type == BC_TYPE_SLIP_X || *bedge_type == BC_TYPE_SLIP_Y || *bedge_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      sP[edge * DG_NPF + i] = s[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      sP[edge * DG_NPF + i] = ps2d_custom_bc_get_ls(*bedge_type, x[fmask_ind], y[fmask_ind], s[fmask_ind]);
    }
  }
}
