inline void zero_npf1(double *npf0) {
  for(int i = 0; i < 3 * DG_NPF; i++) {
    npf0[i] = 0.0;
  }
}
