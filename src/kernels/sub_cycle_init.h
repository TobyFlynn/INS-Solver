inline void sub_cycle_init(const double *u, const double *v, double *u_l, double *v_l) {
  for(int i = 0; i < DG_NP; i++) {
    u_l[i] = u[i];
    v_l[i] = v[i];
  }
}