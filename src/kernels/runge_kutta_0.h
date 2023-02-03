inline void runge_kutta_0(const int *stage, const double *dt, const double *val,
                          const double *rk0, const double *rk1, double *rkQ) {
  if(*stage == -1) {
    for(int i = 0; i < DG_NP; i++) {
      rkQ[i] = val[i];
    }
  } else if(*stage == 0) {
    for(int i = 0; i < DG_NP; i++) {
      rkQ[i] = val[i] + (*dt) * rk0[i];
    }
  } else {
    for(int i = 0; i < DG_NP; i++) {
      rkQ[i] = val[i] + 0.5 * (*dt) * (rk0[i] / 4.0 + rk1[i] / 4.0);
    }
  }
}
