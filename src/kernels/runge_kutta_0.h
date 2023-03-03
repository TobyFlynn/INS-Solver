inline void runge_kutta_0(const int *stage, const DG_FP *dt, const DG_FP *val,
                          const DG_FP *rk0, const DG_FP *rk1, DG_FP *rkQ) {
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
