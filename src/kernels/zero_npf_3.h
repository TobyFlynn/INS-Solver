inline void zero_npf_3(DG_FP *x, DG_FP *y, DG_FP *z) {
  for(int i = 0; i < DG_NUM_FACES * DG_NPF; i++) {
    x[i] = 0.0;
    y[i] = 0.0;
    z[i] = 0.0;
  }
}
