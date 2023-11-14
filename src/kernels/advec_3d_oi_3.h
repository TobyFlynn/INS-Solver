inline void advec_3d_oi_3(const DG_FP *nx, const DG_FP *ny, const DG_FP *nz, const DG_FP *sJ,
                          const DG_FP *geof, const DG_FP *mU, const DG_FP *mV, const DG_FP *mW, 
                          const DG_FP *mVal, DG_FP *pVal) {
  for(int i = 0; i < DG_NUM_FACES * DG_CUB_SURF_3D_NP; i++) {
    const int face_ind = i/DG_CUB_SURF_3D_NP;
    const DG_FP _nx = nx[face_ind];
    const DG_FP _ny = ny[face_ind];
    const DG_FP _nz = nz[face_ind];
    const DG_FP _fscale = sJ[face_ind] / geof[J_IND];

    DG_FP flux0 = _nx * mU[i] + _ny * mV[i] + _nz * mW[i];
    DG_FP flux1 = 0.5 * (mVal[i] + pVal[i]);
    DG_FP flux2 = fabs(flux0);
    DG_FP flux3 = mVal[i] - pVal[i];

    pVal[i] = _fscale * (flux0 * flux1 + 0.5 * flux2 * flux3);
  }
}