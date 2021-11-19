inline void poisson_mf2_faces(const double *uL, const double *opL, double *rhsL,
                              const double *uR, const double *opR, double *rhsR) {
  for(int m = 0; m < DG_NP; m++) {
    int ind = m * DG_NP;
    double valL = 0.0;
    double valR = 0.0;
    for(int n = 0; n < DG_NP; n++) {
      valL += opL[ind + n] * uR[n];
      valR += opR[ind + n] * uL[n];
    }
    rhsL[m] += valL;
    rhsR[m] += valR;
  }
}
