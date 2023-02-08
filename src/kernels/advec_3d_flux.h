inline void advec_3d_flux(const int *faceNum, const int *fmaskL_corrected,
                          const int *fmaskR_corrected, const double *nx,
                          const double *ny, const double *nz,
                          const double *fscale, const double **val,
                          const double **u, const double **v, const double **w,
                          double **flux) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];

  for(int i = 0; i < DG_NPF; i++) {
    int find = faceNum[0] * DG_NPF + i;
    double flux0 = nx[0] * u[0][fmaskL[i]] + ny[0] * v[0][fmaskL[i]] + nz[0] * w[0][fmaskL[i]];
    double flux1 = 0.5 * (val[0][fmaskL[i]] + val[1][fmaskR_corrected[i]]);
    double flux2 = fabs(flux0);
    double flux3 = val[0][fmaskL[i]] - val[1][fmaskR_corrected[i]];
    double flux4 = flux0 * flux1 + 0.5 * flux2 * flux3;

    flux[0][find] += fscale[0] * flux4;
  }

  for(int i = 0; i < DG_NPF; i++) {
    int find = faceNum[1] * DG_NPF + i;
    double flux0 = nx[1] * u[1][fmaskR[i]] + ny[1] * v[1][fmaskR[i]] + nz[1]* w[1][fmaskR[i]];
    double flux1 = 0.5 * (val[1][fmaskR[i]] + val[0][fmaskL_corrected[i]]);
    double flux2 = fabs(flux0);
    double flux3 = val[1][fmaskR[i]] - val[0][fmaskL_corrected[i]];
    double flux4 = flux0 * flux1 + 0.5 * flux2 * flux3;

    flux[1][find] += fscale[1] * flux4;
  }
}
