inline void viscosity_bc(const int *bedge_type, const int *bedgeNum,
                         const double *t, const double *bc_time, const double *x,
                         const double *y, const double *nx, const double *ny,
                         double *exQ0, double *exQ1) {
  int exInd = *bedgeNum * DG_GF_NP;

  const double PI = 3.141592653589793238463;

  if(*bedge_type == 0 || *bedge_type == 3) {
    // Inflow - BC function dependant on time
    if(*t < *bc_time) {
      for(int i = 0; i < DG_GF_NP; i++) {
        double y1 = y[exInd + i];
        exQ0[exInd + i] += sin(PI * (*t) / (*bc_time * 2.0));
      }
    } else {
      for(int i = 0; i < DG_GF_NP; i++) {
        double y1 = y[exInd + i];
        exQ0[exInd + i] += 1.0;
      }
    }
  } else if(*bedge_type == 1) {
    // Outflow - Natural boundary condition
  } else {
    // Wall - No slip
  }
}
