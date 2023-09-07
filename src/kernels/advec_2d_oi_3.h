inline void advec_2d_oi_3(const DG_FP *nx, const DG_FP *ny, const DG_FP *sJ,
                          const DG_FP *geof, const DG_FP *uM, const DG_FP *vM,
                          const DG_FP *valM, DG_FP *valP) {
  for(int i = 0; i < DG_NUM_FACES * DG_CUB_SURF_2D_NP; i++) {
    const int face_ind = i/DG_CUB_SURF_2D_NP;
    const DG_FP _nx = nx[face_ind];
    const DG_FP _ny = ny[face_ind];
    const DG_FP _fscale = sJ[face_ind] / geof[J_IND];

    DG_FP flux0 = _nx * uM[i] + _ny * vM[i];
    DG_FP flux1 = 0.5 * (valM[i] + valP[i]);
    DG_FP flux2 = fabs(flux0);
    DG_FP flux3 = valM[i] - valP[i];

    valP[i] = _fscale * (flux0 * flux1 + 0.5 * flux2 * flux3);
  }
}