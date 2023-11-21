inline void ls_advec_2d_oi_bc(const int *bedge_type, const int *bedgeNum,
                              const DG_FP *x, const DG_FP *y, const DG_FP *u,
                              const DG_FP *v, const DG_FP *val, DG_FP *uM,
                              DG_FP *vM, DG_FP *valM, DG_FP *valP) {
  const int edge = bedgeNum[0];
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];

  for(int i = 0; i < DG_NPF; i++) {
    const int find = fmask[i];
    uM[edge * DG_NPF + i]   = u[find];
    vM[edge * DG_NPF + i]   = v[find];
    valM[edge * DG_NPF + i] = val[find];
  }

  if(*bedge_type == BC_TYPE_NO_SLIP || *bedge_type == BC_TYPE_SLIP_X || *bedge_type == BC_TYPE_SLIP_Y || *bedge_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      valP[edge * DG_NPF + i] = val[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      valP[edge * DG_NPF + i] = ps2d_custom_bc_get_ls(*bedge_type, x[fmask_ind], y[fmask_ind], val[fmask_ind]);
    }
  }
}
