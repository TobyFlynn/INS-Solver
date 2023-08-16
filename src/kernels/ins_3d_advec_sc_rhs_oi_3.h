inline void ins_3d_advec_sc_rhs_oi_3(const DG_FP *nx, const DG_FP *ny, const DG_FP *nz, const DG_FP *sJ,
                              const DG_FP *geof, const DG_FP *mUs, const DG_FP *mVs, const DG_FP *mWs, 
                              const DG_FP *pUs, const DG_FP *pVs, const DG_FP *pWs, const DG_FP *mUb, 
                              const DG_FP *mVb, const DG_FP *mWb, DG_FP *pUb, DG_FP *pVb, DG_FP *pWb) {
  for(int i = 0; i < DG_NUM_FACES * DG_CUB_SURF_3D_NP; i++) {
    const int face_ind = i/DG_CUB_SURF_3D_NP;
    const DG_FP _nx = nx[face_ind];
    const DG_FP _ny = ny[face_ind];
    const DG_FP _nz = nz[face_ind];
    const DG_FP _fscale = sJ[face_ind] / geof[J_IND];

    const DG_FP _mUs = mUs[i];
    const DG_FP _mVs = mVs[i];
    const DG_FP _mWs = mWs[i];
    const DG_FP _pUs = pUs[i];
    const DG_FP _pVs = pVs[i];
    const DG_FP _pWs = pWs[i];
    const DG_FP _mUb = mUb[i];
    const DG_FP _mVb = mVb[i];
    const DG_FP _mWb = mWb[i];
    const DG_FP _pUb = pUb[i];
    const DG_FP _pVb = pVb[i];
    const DG_FP _pWb = pWb[i];

    // Check whether it shoudl be s instead of b here (or the max of them)
    const DG_FP velL = _nx * _mUb + _ny * _mVb + _nz * _mWb;
    const DG_FP velR = _nx * _pUb + _ny * _pVb + _nz * _pWb;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    pUb[i] = 0.5 * _fscale * (_nx * (_mUb * _mUs + _pUb * _pUs) + _ny * (_mUb * _mVs + _pUb * _pVs)
          + _nz * (_mUb * _mWs + _pUb * _pWs) + maxvel * (_mUs - _pUs));
    pVb[i] = 0.5 * _fscale * (_nx * (_mVb * _mUs + _pVb * _pUs) + _ny * (_mVb * _mVs + _pVb * _pVs)
          + _nz * (_mVb * _mWs + _pVb * _pWs) + maxvel * (_mVs - _pVs));
    pWb[i] = 0.5 * _fscale * (_nx * (_mWb * _mUs + _pWb * _pUs) + _ny * (_mWb * _mVs + _pWb * _pVs)
          + _nz * (_mWb * _mWs + _pWb * _pWs) + maxvel * (_mWs - _pWs));
  }
}