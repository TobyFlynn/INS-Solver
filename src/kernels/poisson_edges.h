inline void poisson_edges(const double *uL, const double *opL, double *rhsL,
                          const double *uR, const double *opR, double *rhsR) {
  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind = m + n * DG_NP;
      rhsL[m] += opL[ind] * uR[n];
      rhsR[m] += opR[ind] * uL[n];
    }
  }
}
