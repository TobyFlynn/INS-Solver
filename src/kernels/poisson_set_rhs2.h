inline void poisson_set_rhs2(const double *setBC, double *rhs, double *tau, double *liftx, double *lifty) {
  for(int i = 0; i < 15; i++) {
    // rhs[i] = rhs[i] + (setBC[i] - tau[i]);
    rhs[i] = rhs[i] - setBC[i];
    // rhs[i] = rhs[i] - (liftx[i] - lifty[i]);
    liftx[i] = 0.0;
    lifty[i] = 0.0;
    tau[i] = 0.0;
  }
}
