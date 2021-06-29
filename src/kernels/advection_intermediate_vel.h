inline void advection_intermediate_vel(const double *a0, const double *a1,
                                       const double *b0, const double *b1,
                                       const double *g0, const double *dt,
                                       const double *q0, const double *q1,
                                       const double *q0Old, const double *q1Old,
                                       const double *N0, const double *N1,
                                       const double *N0Old, const double *N1Old,
                                       double *q0T, double *q1T) {
  for(int i = 0; i < 15; i++) {
    // q0T[i] = ((*a0 * q0[i] + *a1 * q0Old[i]) - *dt * (*b0 * N0[i] + *b1 * N0Old[i])) / *g0;
    // q1T[i] = ((*a0 * q1[i] + *a1 * q1Old[i]) - *dt * (*b0 * N1[i] + *b1 * N1Old[i])) / *g0;
    // q0T[i] = ((*a0 * q0[i] + *a1 * q0Old[i]) - *dt * (*b0 * N0[i] + *b1 * N0Old[i]));
    // q1T[i] = ((*a0 * q1[i] + *a1 * q1Old[i]) - *dt * (*b0 * N1[i] + *b1 * N1Old[i]));
    q0T[i] = *a0 * q0[i] + *a1 * q0Old[i] + *dt * (*b0 * N0[i] + *b1 * N0Old[i]);
    q1T[i] = *a0 * q1[i] + *a1 * q1Old[i] + *dt * (*b0 * N1[i] + *b1 * N1Old[i]);
  }
}
