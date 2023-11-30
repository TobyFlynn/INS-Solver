inline void diff_3d_1(const int *faceNum, const DG_FP *nx, const DG_FP *ny, 
                      const DG_FP *nz, const DG_FP *fscale, const DG_FP *val, 
                      DG_FP *flux_x, DG_FP *flux_y, DG_FP *flux_z) {
  // Work out which edge for each element
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;

  const DG_FP _nx = *nx;
  const DG_FP _ny = *ny;
  const DG_FP _nz = *nz;
  const DG_FP _fscale = *fscale;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    const DG_FP _val = val[fmask_ind];
    flux_x[fInd + i] += _fscale * _nx * _val;
    flux_y[fInd + i] += _fscale * _ny * _val;
    flux_z[fInd + i] += _fscale * _nz * _val;
  }
}
