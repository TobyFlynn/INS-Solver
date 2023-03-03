inline void ins_pressure_bc_2d(const DG_FP *t, const int *bedge_type,
                               const int *bedgeNum, const DG_FP *x,
                               const DG_FP *y, const DG_FP *nx,
                               const DG_FP *ny, const DG_FP *N0,
                               const DG_FP *N1, const DG_FP *gradCurlVel0,
                               const DG_FP *gradCurlVel1, DG_FP *dPdN) {
  // Get constants for this element's order
  const int exInd = *bedgeNum * DG_GF_NP;
  const DG_FP PI = 3.141592653589793238463;

  if(*bedge_type == 0 || *bedge_type == 2) {
    // Inflow or Wall
    for(int i = 0; i < DG_GF_NP; i++) {
      int ind = exInd + i;
      DG_FP res1 = -N0[ind] - gradCurlVel1[ind] / r_ynolds;
      DG_FP res2 = -N1[ind] + gradCurlVel0[ind] / r_ynolds;
      dPdN[ind] += nx[ind] * res1 + ny[ind] * res2;
    }
  }
}
