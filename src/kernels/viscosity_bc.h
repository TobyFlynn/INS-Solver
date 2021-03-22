// TODO double check these
inline void viscosity_bc(const int *bedge_type, const int *bedgeNum,
                         const double *t, const double *x, const double *y,
                         double *exQ0, double *exQ1) {
  int exInd = 0;
  if(*bedgeNum == 1) {
    exInd = 7;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 7;
  }

  if(*bedge_type == 0) {
    // Inflow - BC function dependant on time
    const double PI = 3.141592653589793238463;
    for(int i = 0; i < 7; i++) {
      double y1 = y[exInd + i];
      exQ0[exInd + i] += pow(0.41, -2.0) * sin((PI * (*t)) / 8.0) * 6.0 * y1 * (0.41 - y1);
      // exQ0[exInd + i] += pow(0.41, -2.0) * cos((PI * *t) / 8.0);
      // exQ0[exInd + i] += 0.1;
      // exQ1[exInd + i] += bc_v;
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
    // for(int i = 0; i < 7; i++) {
    //   int qInd = fmask[i];
    //   exQ0[exInd + i] += q0[qInd] - 2 * (nx[exInd + i] * q0[qInd] + ny[exInd + i] * q1[qInd]) * nx[exInd + i];
    //   exQ1[exInd + i] += q1[qInd] - 2 * (nx[exInd + i] * q0[qInd] + ny[exInd + i] * q1[qInd]) * ny[exInd + i];
    // }
  }
}
