inline void viscosity_bc(const int *bedge_type, const int *bedgeNum,
                         const double *t, const double *x, const double *y,
                         const double *nx, const double *ny, double *exQ0,
                         double *exQ1) {
  int exInd = *bedgeNum * DG_GF_NP;

  const double PI = 3.141592653589793238463;

  if(*bedge_type == 0) {
    // Inflow - BC function dependant on time
    for(int i = 0; i < DG_GF_NP; i++) {
      double y1 = y[exInd + i];
      exQ0[exInd + i] += (3.0 * PI / 8.0) * sin((PI * (*t)) / 8.0) * y1 * (1.0 - y1);
    }
  } else if(*bedge_type == 1) {
    // Outflow - Natural boundary condition
    // for(int i = 0; i < 7; i++) {
    //   int qInd = fmask[i];
    //   exQ0[exInd + i] += q0[qInd];
    //   exQ1[exInd + i] += q1[qInd];
    // }
  } else {
    // Wall - No slip
  }
}
