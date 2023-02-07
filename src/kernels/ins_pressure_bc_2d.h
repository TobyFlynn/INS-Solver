inline void ins_pressure_bc_2d(const double *t, const int *bedge_type,
                               const int *bedgeNum, const double *x,
                               const double *y, const double *nx,
                               const double *ny, const double *N0,
                               const double *N1, const double *gradCurlVel0,
                               const double *gradCurlVel1, const double *rho,
                               double *dPdN) {
  // Get constants for this element's order
  const int exInd = *bedgeNum * DG_GF_NP;
  const double PI = 3.141592653589793238463;

  if(*bedge_type == 0 || *bedge_type == 2) {
    // Inflow or Wall
    for(int i = 0; i < DG_GF_NP; i++) {
      int ind = exInd + i;
      double res1 = -N0[ind] - gradCurlVel1[ind] / (r_ynolds * rho[ind]);
      double res2 = -N1[ind] + gradCurlVel0[ind] / (r_ynolds * rho[ind]);
      dPdN[ind] += nx[ind] * res1 + ny[ind] * res2;
    }
  }
}
