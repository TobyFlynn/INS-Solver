inline void zero_npf_3_sp(float *x, float *y, float *z) {
  for(int i = 0; i < DG_NUM_FACES * DG_NPF; i++) {
    x[i] = 0.0f;
    y[i] = 0.0f;
    z[i] = 0.0f;
  }
}
