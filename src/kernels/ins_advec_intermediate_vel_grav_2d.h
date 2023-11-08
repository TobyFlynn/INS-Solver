inline void ins_advec_intermediate_vel_grav_2d(const DG_FP *a0, const DG_FP *a1,
                              const DG_FP *b0, const DG_FP *b1, const DG_FP *dt, 
                              const DG_FP *q0, const DG_FP *q1, const DG_FP *q0Old,
                              const DG_FP *q1Old, const DG_FP *N0, const DG_FP *N1, 
                              const DG_FP *N0Old, const DG_FP *N1Old, DG_FP *q0T,
                              DG_FP *q1T) {
  const DG_FP grav_term = 1.0 / (froude * froude);
  for(int i = 0; i < DG_NP; i++) {
    q0T[i] = (*a0 * q0[i] + *a1 * q0Old[i]) - *dt * (*b0 * N0[i] + *b1 * N0Old[i]);
    q1T[i] = (*a0 * q1[i] + *a1 * q1Old[i]) - *dt * (*b0 * (N1[i] + grav_term) + *b1 * (N1Old[i] + grav_term));
  }
}
