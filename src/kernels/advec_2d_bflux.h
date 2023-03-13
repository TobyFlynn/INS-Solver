inline void advec_2d_bflux(const int *p, const int *edgeNum, const DG_FP *nx,
                           const DG_FP *ny, const DG_FP *sJ, const DG_FP *val,
                           const DG_FP *u, const DG_FP *v, DG_FP *flux) {
  // Get constants
  const DG_FP *gaussW = &gaussW_g[(*p - 1) * DG_GF_NP];
  const int exInd = *edgeNum * DG_GF_NP;

  for(int i = 0; i < DG_GF_NP; i++) {
    int ind = exInd + i;

    DG_FP flux0 = nx[ind] * u[ind] + ny[ind] * v[ind];
    DG_FP flux1 = val[ind];
    DG_FP flux2 = fabs(flux0);
    DG_FP flux3 = 0.0;

    flux[ind] += gaussW[i] * sJ[ind] * (flux0 * flux1 + 0.5 * flux2 * flux3);
  }
}
