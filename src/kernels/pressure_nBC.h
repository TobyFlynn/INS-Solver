inline void pressure_nBC(const double *nx, const double *ny, const double *dPdN,
                         double *nBCx, double *nBCy,
                         const double *N0, const double *N1,
                         const double *gradCurlVel0, const double *gradCurlVel1) {
  for(int i = 0; i < 15; i++) {
    // nBCx[i] = nx[i] * dPdN[i];
    // nBCy[i] = ny[i] * dPdN[i];
    int fInd = FMASK[i];
    double res1 = N0[fInd] + nu * gradCurlVel1[fInd];
    double res2 = N1[fInd] - nu * gradCurlVel0[fInd];
    nBCx[i] = -nx[i] * res1;
    nBCy[i] = -ny[i] * res2;
  }
}
