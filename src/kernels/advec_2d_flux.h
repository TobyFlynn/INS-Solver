inline void advec_2d_flux(const int **p, const int *edgeNum, const bool *rev,
                          const DG_FP **nx, const DG_FP **ny,
                          const DG_FP **sJ, const DG_FP **val,
                          const DG_FP **u, const DG_FP **v, DG_FP **flux) {
  // Work out which edge for each element
  const int edgeL = edgeNum[0];
  const int edgeR = edgeNum[1];
  const bool reverse = *rev;
  // Get constants
  const DG_FP *gaussW = &gaussW_g[(p[0][0] - 1) * DG_GF_NP];
  const int exIndL = edgeL * DG_GF_NP;
  const int exIndR = edgeR * DG_GF_NP;

  for(int i = 0; i < DG_GF_NP; i++) {
    int lInd = exIndL + i;
    int rInd;
    int rWInd;
    if(reverse) {
      rInd = exIndR + DG_GF_NP - 1 - i;
      rWInd = DG_GF_NP - 1 - i;
    } else {
      rInd = exIndR + i;
      rWInd = i;
    }
    DG_FP flux0L = nx[0][lInd] * u[0][lInd] + ny[0][lInd] * v[0][lInd];
    DG_FP flux1L = 0.5 * (val[0][lInd] + val[1][rInd]);
    DG_FP flux2L = fabs(flux0L);
    DG_FP flux3L = val[0][lInd] - val[1][rInd];

    DG_FP flux0R = nx[1][rInd] * u[1][rInd] + ny[1][rInd] * v[1][rInd];
    DG_FP flux1R = 0.5 * (val[1][rInd] + val[0][lInd]);
    DG_FP flux2R = fabs(flux0R);
    DG_FP flux3R = val[1][rInd] - val[0][lInd];

    flux[0][lInd] += gaussW[i] * sJ[0][lInd] * (flux0L * flux1L + 0.5 * flux2L * flux3L);
    flux[1][rInd] += gaussW[rWInd] * sJ[1][rInd] * (flux0R * flux1R + 0.5 * flux2R * flux3R);
  }
}
