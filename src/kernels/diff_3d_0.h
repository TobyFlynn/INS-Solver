inline void diff_3d_0(const int *faceNum, const int *fmaskL_corrected,
                      const int *fmaskR_corrected, const DG_FP *nx,
                      const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                      const DG_FP **val, DG_FP **flux_x, DG_FP **flux_y,
                      DG_FP **flux_z) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];
  const int fIndL = faceNum[0] * DG_NPF;
  const int fIndR = faceNum[1] * DG_NPF;

  const DG_FP _nxL = nx[0];
  const DG_FP _nyL = ny[0];
  const DG_FP _nzL = nz[0];
  const DG_FP _fscaleL = fscale[0];
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];
    const DG_FP avg_val = 0.5 * (val[0][fmaskL_ind] + val[1][fmaskR_ind]);

    flux_x[0][fIndL + i] = _fscaleL * _nxL * avg_val;
    flux_y[0][fIndL + i] = _fscaleL * _nyL * avg_val;
    flux_z[0][fIndL + i] = _fscaleL * _nzL * avg_val;
  }

  const DG_FP _nxR = nx[1];
  const DG_FP _nyR = ny[1];
  const DG_FP _nzR = nz[1];
  const DG_FP _fscaleR = fscale[1];
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL_corrected[i];
    const int fmaskR_ind = fmaskR[i];
    const DG_FP avg_val = 0.5 * (val[0][fmaskL_ind] + val[1][fmaskR_ind]);

    flux_x[1][fIndR + i] = _fscaleR * _nxR * avg_val;
    flux_y[1][fIndR + i] = _fscaleR * _nyR * avg_val;
    flux_z[1][fIndR + i] = _fscaleR * _nzR * avg_val;
  }
}
