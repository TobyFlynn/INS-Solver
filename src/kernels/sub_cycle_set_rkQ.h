inline void sub_cycle_set_rkQ(const int *stage, const double *dt, 
                              const double *u, const double *v,
                              const double *rk0_u, const double *rk0_v,
                              const double *rk1_u, const double *rk1_v,
                              double *rkQ_u, double *rkQ_v) {
  if(*stage == -1) {
    for(int i = 0; i < DG_NP; i++) {
      rkQ_u[i] = u[i];
      rkQ_v[i] = v[i];
    }
  } else if(*stage == 0) {
    for(int i = 0; i < DG_NP; i++) {
      rkQ_u[i] = u[i] + (*dt) * rk0_u[i];
      rkQ_v[i] = v[i] + (*dt) * rk0_v[i];
    }
  } else {
    for(int i = 0; i < DG_NP; i++) {
      rkQ_u[i] = u[i] + 0.5 * (*dt) * (rk0_u[i] / 4.0 + rk1_u[i] / 4.0);
      rkQ_v[i] = v[i] + 0.5 * (*dt) * (rk0_v[i] / 4.0 + rk1_v[i] / 4.0);
    }
  }
}
