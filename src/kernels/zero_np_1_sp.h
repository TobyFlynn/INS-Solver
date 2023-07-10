inline void zero_np_1_sp(float *x) {
  for(int i = 0; i < DG_NP; i++) {
    x[i] = 0.0f;
  }
}
