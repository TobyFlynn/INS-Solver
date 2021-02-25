inline void pressure_bc(const int *bedge_type, const int *bedgeNum, const double *t,
                        const double *x, const double *y,
                        const double *nx, const double *ny,
                        const double *N0, const double *N1,
                        const double *gradCurlVel0, const double *gradCurlVel1,
                        double *dPdN) {
  int exInd = 0;
  if(*bedgeNum == 1) {
    exInd = 5;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 5;
  }

  int *fmask;

  if(*bedgeNum == 0) {
    fmask = FMASK;
  } else if(*bedgeNum == 1) {
    fmask = &FMASK[5];
  } else {
    fmask = &FMASK[2 * 5];
  }

  if(*bedge_type == 0 || *bedge_type == 2) {
    // Inflow or Wall
    for(int i = 0; i < 5; i++) {
      int fInd = fmask[i];
      double res1 = N0[fInd] + nu * gradCurlVel1[fInd];
      double res2 = N1[fInd] - nu * gradCurlVel0[fInd];
      dPdN[exInd + i] += nx[exInd + i] * res1 + ny[exInd + i] * res2;
    }
  }

  if(*bedge_type == 0) {
    // Inflow
    const double PI = 3.141592653589793238463;
    for(int i = 0; i < 5; i++) {
      double y1 = y[fmask[i]];
      double bcdUndt = -pow(0.41, -2.0) * (PI/8.0) * cos((PI * *t) / 8.0) * 6.0 * y1 * (0.41 - y1);
      // printf("%g\n", bcdUndt);
      // double bcdUndt = -pow(0.41, -2.0) * (PI/8.0) * cos((PI * *t) / 8.0);
      dPdN[exInd + i] += bcdUndt;
    }
  }
}
