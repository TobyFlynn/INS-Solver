inline void advection_intermediate_vel(const double *a0, const double *a1,
                                       const double *b0, const double *b1,
                                       const double *g0, const double *dt,
                                       const double *q0, const double *q1,
                                       const double *q0Old, const double *q1Old,
                                       const double *N0, const double *N1,
                                       const double *N0Old, const double *N1Old,
                                       const double *st_x, const double *st_y,
                                       const double *st_xOld, const double *st_yOld,
                                       double *q0T, double *q1T) {
  double gravity = 1.0 / (froude * froude);
  for(int i = 0; i < DG_NP; i++) {
    q0T[i] = *a0 * q0[i] + *a1 * q0Old[i] + *dt * (*b0 * N0[i] + *b1 * N0Old[i]);
    q1T[i] = *a0 * q1[i] + *a1 * q1Old[i] + *dt * (*b0 * N1[i] + *b1 * N1Old[i]);
    // q1T[i] = *a0 * q1[i] + *a1 * q1Old[i] + *dt * (*b0 * (N1[i] - 1.0 / (froude * froude)) + *b1 * (N1Old[i] - 1.0 / (froude * froude)));

    // q0T[i] = *a0 * q0[i] + *a1 * q0Old[i] + *dt * (*b0 * (N0[i]) + *b1 * (N0Old[i]));
    // q1T[i] = *a0 * q1[i] + *a1 * q1Old[i] + *dt * (*b0 * (N1[i] - gravity) + *b1 * (N1Old[i] - gravity));

    // q0T[i] = *a0 * q0[i] + *a1 * q0Old[i] + *dt * (*b0 * (N0[i] + st_x[i]) + *b1 * (N0Old[i] + st_xOld[i]));
    // q1T[i] = *a0 * q1[i] + *a1 * q1Old[i] + *dt * (*b0 * (N1[i] - gravity + st_y[i]) + *b1 * (N1Old[i] - gravity + st_yOld[i]));
  }
}
