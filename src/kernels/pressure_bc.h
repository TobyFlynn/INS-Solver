inline void pressure_bc(const int *bedge_type, const int *bedgeNum,
                        const double *t, const int *problem, const double *x,
                        const double *y, const double *nx, const double *ny,
                        const double *N0, const double *N1,
                        const double *gradCurlVel0, const double *gradCurlVel1,
                        double *dPdN) {
  int exInd = *bedgeNum * DG_NPF;
  int *fmask = &FMASK[*bedgeNum * DG_NPF];

  const double PI = 3.141592653589793238463;

  if(*problem == 0) {
    if(*bedge_type == 0 || *bedge_type == 2 || *bedge_type == 3) {
      // Inflow or Wall
      for(int i = 0; i < DG_NPF; i++) {
        int fInd = fmask[i];
        double res1 = -N0[fInd] - nu * gradCurlVel1[fInd];
        double res2 = -N1[fInd] + nu * gradCurlVel0[fInd];
        dPdN[exInd + i] += nx[exInd + i] * res1 + ny[exInd + i] * res2;
      }
    }

    if(*bedge_type == 0) {
      // Inflow
      for(int i = 0; i < DG_NPF; i++) {
        double y1 = y[fmask[i]];
        double bcdUndt = -pow(0.41, -2.0) * (PI/8.0) * cos((PI * *t) / 8.0) * 6.0 * y1 * (0.41 - y1);
        dPdN[exInd + i] -= bcdUndt;
      }
    }
  } else {
    if(*bedge_type == 0) {
      // Inflow
      for(int i = 0; i < DG_NPF; i++) {
        int fInd = fmask[i];
        double res1 = -N0[fInd] - nu * gradCurlVel1[fInd];
        double res2 = -N1[fInd] + nu * gradCurlVel0[fInd];
        dPdN[exInd + i] += nx[exInd + i] * res1 + ny[exInd + i] * res2;

        double y1 = y[fmask[i]];
        double x1 = x[fmask[i]];
        double nx1 = nx[exInd + i];
        double ny1 = ny[exInd + i];
        double bcdUndt = -nu * 4.0 * PI * PI * (-nx1 * sin(2.0 * PI * y1) + ny1 * sin(2.0 * PI * x1))
                          * exp(-nu * 4.0 * PI * PI * *t);
        dPdN[exInd + i] -= bcdUndt;
      }
    }
  }
}
