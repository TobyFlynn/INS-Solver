inline void zero_npf_2(DG_FP *npf0, DG_FP *npf1) {
  for(int i = 0; i < DG_NUM_FACES * DG_NPF; i++) {
    npf0[i] = 0.0;
    npf1[i] = 0.0;
  }
}
