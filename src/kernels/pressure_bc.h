inline void pressure_bc(const int *bedge_type, const int *bedgeNum,
                        const int *p, const double *t, const double *x,
                        const double *y, const double *nx, const double *ny,
                        const double *N0, const double *N1,
                        const double *gradCurlVel0, const double *gradCurlVel1,
                        const double *rho, double *dPdN) {
  // Get constants for this element's order
  const int dg_npf = DG_CONSTANTS[(*p - 1) * 5 + 1];
  const int *fmask = &FMASK[(*p - 1) * 3 * DG_NPF];
  fmask = &fmask[*bedgeNum * dg_npf];
  const int exInd = *bedgeNum * dg_npf;

  const double PI = 3.141592653589793238463;

  if(*bedge_type == 0 || *bedge_type == 2 || *bedge_type == 3) {
    // Inflow or Wall
    for(int i = 0; i < dg_npf; i++) {
      int fInd = fmask[i];
      double res1 = -N0[fInd] - gradCurlVel1[fInd] / (reynolds * rho[fInd]);
      double res2 = -N1[fInd] + gradCurlVel0[fInd] / (reynolds * rho[fInd]);
      dPdN[exInd + i] += nx[exInd + i] * res1 + ny[exInd + i] * res2;
    }
  }

  if(*bedge_type == 0) {
    // Inflow
    for(int i = 0; i < dg_npf; i++) {
      double y1 = y[fmask[i]];
      double bcdUndt = -(PI/8.0) * (3.0 * PI / 8.0) * cos((PI * *t) / 8.0) * y1 * (1.0 - y1);
      dPdN[exInd + i] -= bcdUndt;
    }
  }
}
