inline void poisson_edges(const double *uL, const double *opL, double *rhsL,
                          const double *uR, const double *opR, double *rhsR) {
  for(int m = 0; m < 3; m++) {
    int ind = m * 3;
    for(int n = 0; n < 3; n++) {
      rhsL[m] += opL[ind + n] * uR[n];
      rhsR[m] += opR[ind + n] * uL[n];
    }
  }
}
