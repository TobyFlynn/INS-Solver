inline void ins_vis_bc_2d_dam(const DG_FP *t, const DG_FP *g0, const int *bedge_type,
                            const int *bedgeNum, const DG_FP *x, const DG_FP *y,
                            const DG_FP *nx, const DG_FP *ny, const DG_FP *velTT0,
                            const DG_FP *velTT1, DG_FP *exQ0, DG_FP *exQ1) {
  int exInd = *bedgeNum * DG_GF_NP;
  const DG_FP PI = 3.141592653589793238463;

  if(*bedge_type == 1) {
    // Outflow - Natural boundary condition
  } else {
    // Wall - slip
    for(int i = 0; i < DG_GF_NP; i++) {
      const DG_FP dot = nx[exInd + i] * velTT0[exInd + i] + ny[exInd + i] * velTT1[exInd + i];
      exQ0[exInd + i] = (velTT0[exInd + i] - dot * nx[exInd + i]) / *g0;
      exQ1[exInd + i] = (velTT1[exInd + i] - dot * ny[exInd + i]) / *g0;
    }
  }
}
