inline void poisson_mf2_faces(const double *uL, const double *opL, double *rhsL,
                              const double *uR, const double *opR, double *rhsR) {
  for(int m = 0; m < DG_NP; m++) {
    int ind = m * DG_NP;
    double val = 0.0;
    for(int n = 0; n < DG_NP; n++) {
      val += opL[ind + n] * uR[n];
    }
    rhsL[m] += val;
  }

  for(int m = 0; m < DG_NP; m++) {
    int ind = m * DG_NP;
    double val = 0.0;
    for(int n = 0; n < DG_NP; n++) {
      val += opR[ind + n] * uL[n];
    }
    rhsR[m] += val;
  }
}
