inline void diff_3d_4(const int *faceNum, const DG_FP *nx, const DG_FP *ny, 
                      const DG_FP *nz, const DG_FP *fscale, const DG_FP *val_x, 
                      const DG_FP *val_y, const DG_FP *val_z, const DG_FP *val, 
                      const DG_FP *vis, DG_FP *flux) {
  // Work out which edge for each element
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;

  const DG_FP _nx = *nx;
  const DG_FP _ny = *ny;
  const DG_FP _nz = *nz;
  const DG_FP _fscale = *fscale;
  const DG_FP pen = DG_ORDER * DG_ORDER * _fscale;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];

    // Zero gradient on boundary
    const DG_FP avg_val_x = 0.5 * val_x[fmask_ind];
    const DG_FP avg_val_y = 0.5 * val_y[fmask_ind];
    const DG_FP avg_val_z = 0.5 * val_z[fmask_ind];
    const DG_FP diff_val = 0.0;

    flux[fInd + i] = _fscale * (_nx * avg_val_x + _ny * avg_val_y + _nz * avg_val_z - pen * vis[fmask_ind] * diff_val);
  }
}
