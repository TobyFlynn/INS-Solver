inline void set_ic(double *q0, double *q1, double *q2) {
  for(int i = 0; i < 15; i++) {
    q0[i] = bc_u;
    q1[i] = bc_v;
    q2[i] = bc_p;
  }
}
