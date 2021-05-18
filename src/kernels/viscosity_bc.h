inline void viscosity_bc(const int *bedge_type, const int *bedgeNum,
                         const double *t, const int *problem, const double *x,
                         const double *y, const double *nx, const double *ny,
                         double *exQ0, double *exQ1) {
  int exInd = 0;
  if(*bedgeNum == 1) {
    exInd = 7;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 7;
  }

  const double PI = 3.141592653589793238463;

  if(*problem == 0) {
    if(*bedge_type == 0) {
      // Inflow - BC function dependant on time
      for(int i = 0; i < 7; i++) {
        double y1 = y[exInd + i];
        exQ0[exInd + i] += pow(0.41, -2.0) * sin((PI * (*t)) / 8.0) * 6.0 * y1 * (0.41 - y1);
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
  } else {
    if(*bedge_type == 0) {
      // Inflow - BC function dependant on time
      for(int i = 0; i < 7; i++) {
        double y1 = y[exInd + i];
        double x1 = x[exInd + i];
        exQ0[exInd + i] += -sin(2.0 * PI * y1) * exp(-nu * 4.0 * PI * PI * *t);
        exQ1[exInd + i] += sin(2.0 * PI * x1) * exp(-nu * 4.0 * PI * PI * *t);
      }
    }

    if(*bedge_type == 1) {
      // Outflow - BC function dependant on time
      for(int i = 0; i < 7; i++) {
        double y1  = y[exInd + i];
        double x1  = x[exInd + i];
        double ny1 = ny[exInd + i];
        double nx1 = nx[exInd + i];
        exQ0[exInd + i] += ny1 * 2.0 * PI * (-cos(2.0 * PI * y1)) * exp(-nu * 4.0 * PI * PI * *t);
        exQ1[exInd + i] += nx1 * 2.0 * PI * cos(2.0 * PI * x1) * exp(-nu * 4.0 * PI * PI * *t);
      }
    }
  }
}
