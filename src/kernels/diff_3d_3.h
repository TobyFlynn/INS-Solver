inline void diff_3d_3(const int *faceNum, const int *fmaskL_corrected,
                      const int *fmaskR_corrected, const DG_FP *nx,
                      const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale, 
                      const DG_FP **val_x, const DG_FP **val_y, const DG_FP **val_z, 
                      const DG_FP **val, const DG_FP **vis, DG_FP **flux) {
  // Work out which edge for each element
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];
  const int fIndL = faceNum[0] * DG_NPF;
  const int fIndR = faceNum[1] * DG_NPF;

  const DG_FP pen = DG_ORDER * DG_ORDER * fmax(fscale[0], fscale[1]);

  const DG_FP _nxL = nx[0];
  const DG_FP _nyL = ny[0];
  const DG_FP _nzL = nz[0];
  const DG_FP _fscaleL = fscale[0];
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];

    const DG_FP avg_val_x = 0.5 * (val_x[0][fmaskL_ind] + val_x[1][fmaskR_ind]);
    const DG_FP avg_val_y = 0.5 * (val_y[0][fmaskL_ind] + val_y[1][fmaskR_ind]);
    const DG_FP avg_val_z = 0.5 * (val_z[0][fmaskL_ind] + val_z[1][fmaskR_ind]);
    const DG_FP diff_val = val[0][fmaskL_ind] - val[1][fmaskR_ind];

    flux[0][fIndL + i] = _fscaleL * (_nxL * avg_val_x + _nyL * avg_val_y + _nzL * avg_val_z - pen * vis[0][fmaskL_ind] * diff_val);
  }

  const DG_FP _nxR = nx[1];
  const DG_FP _nyR = ny[1];
  const DG_FP _nzR = nz[1];
  const DG_FP _fscaleR = fscale[1];
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL_corrected[i];
    const int fmaskR_ind = fmaskR[i];

    const DG_FP avg_val_x = 0.5 * (val_x[0][fmaskL_ind] + val_x[1][fmaskR_ind]);
    const DG_FP avg_val_y = 0.5 * (val_y[0][fmaskL_ind] + val_y[1][fmaskR_ind]);
    const DG_FP avg_val_z = 0.5 * (val_z[0][fmaskL_ind] + val_z[1][fmaskR_ind]);
    const DG_FP diff_val = val[1][fmaskR_ind] - val[0][fmaskL_ind];

    flux[1][fIndR + i] = _fscaleR * (_nxR * avg_val_x + _nyR * avg_val_y + _nzR * avg_val_z - pen * vis[1][fmaskR_ind] * diff_val);
  }
}
