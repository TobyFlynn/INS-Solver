inline void ls_advec_3d_oi_bc(const int *faceNum, const int *bc_type, const DG_FP *x,
                           const DG_FP *y, const DG_FP *z, const DG_FP *nx,
                           const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                           const DG_FP *val, const DG_FP *u, const DG_FP *v, const DG_FP *w,
                           DG_FP *mU, DG_FP *mV, DG_FP *mW, DG_FP *mVal, DG_FP *pVal) {
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    mU[fInd + i]   += u[fmask_ind];
    mV[fInd + i]   += v[fmask_ind];
    mW[fInd + i]   += w[fmask_ind];
    mVal[fInd + i] += val[fmask_ind];
  }

  if(*bc_type == BC_TYPE_NO_SLIP || *bc_type == BC_TYPE_SLIP || *bc_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      pVal[fInd + i] += val[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      pVal[fInd + i] += ps3d_custom_bc_get_ls(*bc_type, x[fmask_ind], y[fmask_ind], z[fmask_ind], val[fmask_ind]);
    }
  }
}
