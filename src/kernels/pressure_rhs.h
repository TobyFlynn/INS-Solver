inline void pressure_rhs(const double *b0, const double *b1, const double *g0,
                         const double *dt, const double *J, const double *sJ,
                         const double *dPdN, double *dPdNOld, double *divVelT) {
  for(int i = 0; i < 15; i++) {
    divVelT[i] = J[i] * (-divVelT[i] * *g0 / *dt);
    dPdNOld[i] = sJ[i] * (*b0 * dPdN[i] + *b1 * dPdNOld[i]);
  }
}
