inline void set_ic(double *q0, double *q1, double *q2, double *q3) {
  for(int i = 0; i < 15; i++) {
    q0[i] = bc_r;
    q1[i] = bc_r * bc_u;
    q2[i] = bc_r * bc_v;
    q3[i] = bc_e;
  }
}
