inline void poisson_mf2_mass(const double *u, const double *op, const double *factor,
                             const double *mm, double *rhs) {
  double mFactor = *factor;
  for(int m = 0; m < 15; m++) {
    int ind = m * 15;
    double val = 0.0;
    for(int n = 0; n < 15; n++) {
      // mm is in column major format while op is in row major
      val += (op[ind + n] + mm[n * 15 + m] * mFactor) * u[n];
    }
    rhs[m] = val;
  }
}
