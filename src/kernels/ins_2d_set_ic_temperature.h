inline void ins_2d_set_ic_temperature(const DG_FP *x, const DG_FP *y,
                                      DG_FP *temperature) {
  for(int i = 0; i < DG_NP; i++) {
    temperature[i] = x[i] / 10.0;
  }
}
