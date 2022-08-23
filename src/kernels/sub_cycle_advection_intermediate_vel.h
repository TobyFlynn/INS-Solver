inline void sub_cycle_advection_intermediate_vel(
                const double *a0, const double *a1,
                const double *q0,const double *q1, 
                const double *q0Old, const double *q1Old, 
                double *q0T, double *q1T) {
  for(int i = 0; i < DG_NP; i++) {
    q0T[i] = *a0 * q0[i] + *a1 * q0Old[i];
    q1T[i] = *a0 * q1[i] + *a1 * q1Old[i];
  }
}
