inline void advec_2d_flux(const int *edgeNum, const bool *rev, const DG_FP *nx,
                          const DG_FP *ny, const DG_FP *fscale,
                          const DG_FP **val, const DG_FP **u, const DG_FP **v,
                          DG_FP **flux) {
  // Work out which edge for each element
  const int edgeL = edgeNum[0];
  const int edgeR = edgeNum[1];
  const bool reverse = *rev;
  const int *fmaskL = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeL * DG_NPF];
  const int *fmaskR = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeR * DG_NPF];

  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = reverse ? fmaskR[DG_NPF - i - 1] : fmaskR[i];

    DG_FP flux0L = nx[0] * u[0][fmaskL_ind] + ny[0] * v[0][fmaskL_ind];
    DG_FP flux1L = 0.5 * (val[0][fmaskL_ind] + val[1][fmaskR_ind]);
    DG_FP flux2L = fabs(flux0L);
    DG_FP flux3L = val[0][fmaskL_ind] - val[1][fmaskR_ind];

    DG_FP flux0R = nx[1] * u[1][fmaskR_ind] + ny[1] * v[1][fmaskR_ind];
    DG_FP flux1R = 0.5 * (val[1][fmaskR_ind] + val[0][fmaskL_ind]);
    DG_FP flux2R = fabs(flux0R);
    DG_FP flux3R = val[1][fmaskR_ind] - val[0][fmaskL_ind];

    const int outL_ind = edgeL * DG_NPF + i;
    const int outR_ind = reverse ? edgeR * DG_NPF + DG_NPF - i - 1 : edgeR * DG_NPF + i;

    flux[0][outL_ind] += fscale[0] * (flux0L * flux1L + 0.5 * flux2L * flux3L);
    flux[1][outR_ind] += fscale[1] * (flux0R * flux1R + 0.5 * flux2R * flux3R);
  }
}
