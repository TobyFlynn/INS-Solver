// TODO double check these
inline void advection_bc(const int *bedge_type, const int *bedgeNum,
                         const double *t, const int *problem,
                         const double *x, const double *y, const double *nu,
                         const double *q0, const double *q1, double *exQ0,
                         double *exQ1) {
  int exInd = *bedgeNum * 4;
  int *fmask = &FMASK[*bedgeNum * 4];

  const double PI = 3.141592653589793238463;

  if(*problem == 0) {
    if(*bedge_type == 0) {
      // Inflow - BC function dependant on time
      for(int i = 0; i < 4; i++) {
        int qInd = fmask[i];
        double y1 = y[qInd];
        exQ0[exInd + i] += pow(1.0, -2.0) * sin((PI * *t) / 8.0) * 6.0 * y1 * (1.0 - y1);
      }
    } else if(*bedge_type == 1) {
      // Outflow - Natural boundary condition
      for(int i = 0; i < 4; i++) {
        int qInd = fmask[i];
        exQ0[exInd + i] += q0[qInd];
        exQ1[exInd + i] += q1[qInd];
      }
    } else {
      // Wall - No slip

    }
  } else {
    if(*bedge_type == 0) {
      // Inflow - BC function dependant on time
      for(int i = 0; i < 4; i++) {
        int qInd = fmask[i];
        double y1 = y[qInd];
        double x1 = x[qInd];
        exQ0[exInd + i] += -sin(2.0 * PI * y1) * exp(-nu[qInd] * 4.0 * PI * PI * *t);
        exQ1[exInd + i] += sin(2.0 * PI * x1) * exp(-nu[qInd] * 4.0 * PI * PI * *t);
      }
    }

    if(*bedge_type == 1) {
      // Outflow - BC function dependant on time
      for(int i = 0; i < 4; i++) {
        int qInd = fmask[i];
        // exQ0[exInd + i] += q0[qInd];
        // exQ1[exInd + i] += q1[qInd];
        double y1 = y[qInd];
        double x1 = x[qInd];
        exQ0[exInd + i] += -sin(2.0 * PI * y1) * exp(-nu[qInd] * 4.0 * PI * PI * *t);
        exQ1[exInd + i] += sin(2.0 * PI * x1) * exp(-nu[qInd] * 4.0 * PI * PI * *t);
      }
    }
  }
}
