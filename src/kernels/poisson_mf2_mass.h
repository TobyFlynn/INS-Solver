inline void poisson_mf2_mass(const double *u, const double *op, const double *factor,
                             const double *mm, double *rhs) {
  double mFactor = *factor;
  for(int m = 0; m < DG_NP; m++) {
    int ind = m * DG_NP;
    double val = 0.0;
    for(int n = 0; n < DG_NP; n++) {
      // mm is in column major format while op is in row major
      val += (op[ind + n] + mm[n * DG_NP + m] * mFactor) * u[n];
    }
    rhs[m] = val;
  }
}
