inline void sub_cycle_update_Q(const double *dt, const double *rk0_u, const double *rk0_v,
                               const double *rk1_u, const double *rk1_v,
                               const double *rk2_u, const double *rk2_v,
                               double *u_l, double *v_l) {
  for(int i = 0; i < DG_NP; i++) {
    double add_u = (*dt) * (rk0_u[i]/ 6.0 + rk1_u[i] / 6.0 + 2.0 * rk2_u[i] / 3.0);
    double add_v = (*dt) * (rk0_v[i]/ 6.0 + rk1_v[i] / 6.0 + 2.0 * rk2_v[i] / 3.0);
    u_l[i] = u_l[i] + add_u;
    v_l[i] = v_l[i] + add_v;
  }
}
