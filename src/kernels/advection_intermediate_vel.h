inline void advection_intermediate_vel(const double *a0, const double *a1,
                                       const double *b0, const double *b1,
                                       const double *dt, const double *q0,
                                       const double *q1, const double *q0Old,
                                       const double *q1Old, const double *N0,
                                       const double *N1, const double *N0Old,
                                       const double *N1Old, const double *sf_x,
                                       const double *sf_y, double *bf_x, double *bf_y,
                                       double *q0T, double *q1T) {
  for(int i = 0; i < DG_NP; i++) {
    q0T[i] = (*a0 * q0[i] + *a1 * q0Old[i]) - *dt * (*b0 * (N0[i] + sf_x[i]) + *b1 * (N0Old[i] + bf_x[i]));
    q1T[i] = (*a0 * q1[i] + *a1 * q1Old[i]) - *dt * (*b0 * (N1[i] + sf_y[i]) + *b1 * (N1Old[i] + bf_y[i]));
    bf_x[i] = sf_x[i];
    bf_y[i] = sf_y[i];
  }
}
