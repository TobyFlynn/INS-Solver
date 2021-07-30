inline void set_ic(double *q0, double *q1) {
  for(int i = 0; i < 10; i++) {
    q0[i] = ic_u;
    q1[i] = ic_v;
  }
}
