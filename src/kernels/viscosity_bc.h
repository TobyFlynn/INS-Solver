inline void viscosity_bc(const int *bedge_type, const int *bedgeNum,
                         const double *t, const int *problem, const double *x,
                         const double *y, const double *nx, const double *ny,
                         double *exQ0, double *exQ1) {
  int exInd = *bedgeNum * DG_GF_NP;

  const double PI = 3.141592653589793238463;

  if(*bedge_type == 0 || *bedge_type == 3) {
    // Inflow - BC function dependant on time
    for(int i = 0; i < DG_GF_NP; i++) {
      double y1 = y[exInd + i];
      exQ0[exInd + i] += sin((PI * (*t)) / 8.0);
    }
  } else if(*bedge_type == 1) {
    // Outflow - Natural boundary condition
  } else {
    // Wall - No slip
  }
}
