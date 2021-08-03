inline void zero_npf(double *npf0, double *npf1) {
  for(int i = 0; i < 3 * 4; i++) {
    npf0[i] = 0.0;
    npf1[i] = 0.0;
  }
}
