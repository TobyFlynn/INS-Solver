inline void poisson_mf2_faces(const double *uL, const double *opL, double *rhsL,
                              const double *uR, const double *opR, double *rhsR) {
  for(int m = 0; m < 15; m++) {
    int ind = m * 15;
    double val = 0.0;
    for(int n = 0; n < 15; n++) {
      val += opL[ind + n] * uR[n];
    }
    rhsL[m] += val;
  }

  for(int m = 0; m < 15; m++) {
    int ind = m * 15;
    double val = 0.0;
    for(int n = 0; n < 15; n++) {
      val += opR[ind + n] * uL[n];
    }
    rhsR[m] += val;
  }
}
