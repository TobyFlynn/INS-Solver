inline void p_multigrid_runge_kutta_1(const DG_FP *dt, DG_FP *val, const DG_FP *rk0,
                          const DG_FP *rk1, const DG_FP *rk2, const DG_FP *b) {
  for(int i = 0; i < DG_NP; i++) {
    DG_FP add = (*dt) * ((b[i] - rk0[i])/ 6.0 + (b[i] - rk1[i]) / 6.0 + 2.0 * (b[i] - rk2[i]) / 3.0);
    val[i] += add;
  }
}
