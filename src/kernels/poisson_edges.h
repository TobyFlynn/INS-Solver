inline void poisson_edges(const double *uL, const double *opL, double *rhsL,
                          const double *uR, const double *opR, double *rhsR) {
  for(int m = 0; m < DG_NP; m++) {
    int ind = m * DG_NP;
    for(int n = 0; n < DG_NP; n++) {
      rhsL[m] += opL[ind + n] * uR[n];
      rhsR[m] += opR[ind + n] * uL[n];
    }
  }
}
