inline void advec_3d_flux(const int *faceNum, const int *fmaskL_corrected,
                          const int *fmaskR_corrected, const DG_FP *nx,
                          const DG_FP *ny, const DG_FP *nz,
                          const DG_FP *fscale, const DG_FP **val,
                          const DG_FP **u, const DG_FP **v, const DG_FP **w,
                          DG_FP **flux) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int faceNumL = faceNum[0];
  const int faceNumR = faceNum[1];
  const int *fmaskL = &fmask[faceNumL * DG_NPF];
  const int *fmaskR = &fmask[faceNumR * DG_NPF];

  for(int i = 0; i < DG_NPF; i++) {
    const int find = faceNum[0] * DG_NPF + i;
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];
    DG_FP flux0 = nx[0] * u[0][fmaskL_ind] + ny[0] * v[0][fmaskL_ind] + nz[0] * w[0][fmaskL_ind];
    DG_FP flux1 = 0.5 * (val[0][fmaskL_ind] + val[1][fmaskR_ind]);
    DG_FP flux2 = fabs(flux0);
    DG_FP flux3 = val[0][fmaskL_ind] - val[1][fmaskR_ind];
    DG_FP flux4 = flux0 * flux1 + 0.5 * flux2 * flux3;
    // DG_FP flux4 = flux0 * flux3;

    flux[0][find] += fscale[0] * flux4;
    // flux[0][find] += fscale[0] * (flux0 * val[0][fmaskL[i]] - flux4);
  }

  for(int i = 0; i < DG_NPF; i++) {
    const int find = faceNum[1] * DG_NPF + i;
    const int fmaskL_ind = fmaskL_corrected[i];
    const int fmaskR_ind = fmaskR[i];
    DG_FP flux0 = nx[1] * u[1][fmaskR_ind] + ny[1] * v[1][fmaskR_ind] + nz[1]* w[1][fmaskR_ind];
    DG_FP flux1 = 0.5 * (val[1][fmaskR_ind] + val[0][fmaskL_ind]);
    DG_FP flux2 = fabs(flux0);
    DG_FP flux3 = val[1][fmaskR_ind] - val[0][fmaskL_ind];
    DG_FP flux4 = flux0 * flux1 + 0.5 * flux2 * flux3;
    // DG_FP flux4 = flux0 * flux3;

    flux[1][find] += fscale[1] * flux4;
    // flux[1][find] += fscale[1] * (flux0 * val[1][fmaskR[i]] - flux4);
  }
}
