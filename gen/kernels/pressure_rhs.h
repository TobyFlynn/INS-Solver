inline void pressure_rhs(const double *b0, const double *b1, const double *g0,
                         const double *dt, const double *J, const double *sJ,
                         const double *dPdN, double *dPdNOld, double *divVelT) {
  // double factor = (*g0) / (*dt);
  double factor = 1.0 / (*dt);
  for(int i = 0; i < 10; i++) {
    divVelT[i] = J[i] * (-divVelT[i] * factor);
  }

  for(int i = 0; i < 3 * 4; i++) {
    dPdNOld[i] = sJ[i] * ((*b0) * dPdN[i] + (*b1) * dPdNOld[i]);
  }
}
