inline void set_rkQ(const int *stage, const double *dt, const double *s,
                    const double *rk0, const double *rk1, double *rkQ) {
  if(*stage == -1) {
    for(int i = 0; i < 10; i++) {
      rkQ[i] = s[i];
    }
  } else if(*stage == 0) {
    for(int i = 0; i < 10; i++) {
      rkQ[i] = s[i] + (*dt) * rk0[i];
    }
  } else {
    for(int i = 0; i < 10; i++) {
      rkQ[i] = s[i] + 0.5 * (*dt) * (rk0[i] / 4.0 + rk1[i] / 4.0);
    }
  }
}
