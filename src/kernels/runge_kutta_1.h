inline void runge_kutta_1(const DG_FP *dt, DG_FP *val, const DG_FP *rk0,
                          const DG_FP *rk1, const DG_FP *rk2) {
  for(int i = 0; i < DG_NP; i++) {
    DG_FP add = (*dt) * (rk0[i]/ 6.0 + rk1[i] / 6.0 + 2.0 * rk2[i] / 3.0);
    val[i] = val[i] + add;
  }
}
