inline void ins_3d_advec_oi_3(const DG_FP *nx, const DG_FP *ny, const DG_FP *nz, const DG_FP *sJ,
                              const DG_FP *geof, const DG_FP *mU, const DG_FP *mV, const DG_FP *mW, 
                              DG_FP *pU, DG_FP *pV, DG_FP *pW) {
  for(int i = 0; i < DG_NUM_FACES * DG_CUB_SURF_3D_NP; i++) {
    const int face_ind = i/DG_CUB_SURF_3D_NP;
    const DG_FP _nx = nx[face_ind];
    const DG_FP _ny = ny[face_ind];
    const DG_FP _nz = nz[face_ind];
    const DG_FP _fscale = sJ[face_ind] / geof[J_IND];

    const DG_FP _mU = mU[i];
    const DG_FP _mV = mV[i];
    const DG_FP _mW = mW[i];
    const DG_FP _pU = pU[i];
    const DG_FP _pV = pV[i];
    const DG_FP _pW = pW[i];

    const DG_FP velL = _nx * _mU + _ny * _mV + _nz * _mW;
    const DG_FP velR = _nx * _pU + _ny * _pV + _nz * _pW;
    const DG_FP maxvel = fmax(fabs(velL), fabs(velR));

    pU[i] = 0.5 * _fscale * (_nx * (_mU * _mU + _pU * _pU) + _ny * (_mU * _mV + _pU * _pV)
          + _nz * (_mU * _mW + _pU * _pW) + maxvel * (_mU - _pU));
    pV[i] = 0.5 * _fscale * (_nx * (_mV * _mU + _pV * _pU) + _ny * (_mV * _mV + _pV * _pV)
          + _nz * (_mV * _mW + _pV * _pW) + maxvel * (_mV - _pV));
    pW[i] = 0.5 * _fscale * (_nx * (_mW * _mU + _pW * _pU) + _ny * (_mW * _mV + _pW * _pV)
          + _nz * (_mW * _mW + _pW * _pW) + maxvel * (_mW - _pW));
  }
}