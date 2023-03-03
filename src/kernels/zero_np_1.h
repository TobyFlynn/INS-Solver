inline void zero_np_1(DG_FP *x) {
  for(int i = 0; i < DG_NP; i++) {
    x[i] = 0.0;
  }
}
